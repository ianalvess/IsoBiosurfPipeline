"""Orchestrator for the isolate genomics pipeline (isolados-biosurf).

Currently implements:
  1. config loading (config/samples.yaml);
  2. Docker daemon check;
  3. QC step (fastp), run as a container. Paired-end only.
  4. assembly step (SPAdes --isolate), run as a container.
  5. QUAST: assembly quality metrics (N50, contig count, etc.).
  6. CheckM2: completeness/contamination estimate for the assembled genome.
  7. BioSurfDB search: Prodigal gene prediction on the assembled genome +
     DIAMOND search against a local BioSurfDB database, summarized at three
     levels of the BioSurfDB hierarchy (categories, subclasses, broad
     classes) as CSV tables and bar charts. With --hits, only this report
     is generated from an existing hits TSV.
"""

import sys
from pathlib import Path

import click
import docker
import yaml
from loguru import logger

DOCKER_IMAGES = {
    "qc": "isolados-biosurf/fastp",
    "assembly": "isolados-biosurf/spades",
    "quast": "isolados-biosurf/quast",
    "checkm2": "isolados-biosurf/checkm2",
    "biosurfdb": "isolados-biosurf/biosurfdb",
}


def load_samples(config_path: Path) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)["samples"]


def to_docker_path(path: Path) -> str:
    """Normalize a resolved path to forward slashes (C:/foo/bar) for
    docker-py bind mounts. docker-py does not apply the Windows-path
    translation the docker CLI does, so backslash paths can be
    misinterpreted; forward-slash paths with the drive letter kept work
    correctly against Docker Desktop. No-op path separator normalization
    on non-Windows platforms.
    """
    return path.as_posix()


def check_docker() -> docker.DockerClient | None:
    try:
        client = docker.from_env()
        client.ping()
        return client
    except Exception as exc:  # noqa: BLE001
        logger.error(f"Could not reach the Docker daemon: {exc}")
        return None


def run_qc(
    client: docker.DockerClient,
    sample_id: str,
    r1: Path,
    r2: Path | None,
    output_dir: Path,
) -> dict:
    image = DOCKER_IMAGES["qc"]
    output_dir.mkdir(parents=True, exist_ok=True)
    paired = r2 is not None

    clean_r1_name = f"{sample_id}_R1.clean.fastq.gz"

    container_cmd = [
        "-i", f"/input/{r1.name}",
        "-o", f"/output/{clean_r1_name}",
        "-j", "/output/fastp.json",
        "-h", "/output/fastp.html",
    ]

    clean_reads = {"r1": output_dir / clean_r1_name}

    if paired:
        clean_r2_name = f"{sample_id}_R2.clean.fastq.gz"
        container_cmd += [
            "-I", f"/input/{r2.name}",
            "-O", f"/output/{clean_r2_name}",
        ]
        clean_reads["r2"] = output_dir / clean_r2_name

    logger.info(
        f"Running QC step (fastp) with image '{image}' for sample '{sample_id}' "
        f"({'paired-end' if paired else 'single-end'})..."
    )

    volumes = {
        to_docker_path(r1.parent.resolve()): {"bind": "/input", "mode": "ro"},
        to_docker_path(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
    }

    logs = client.containers.run(
        image,
        command=container_cmd,
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success(f"QC step finished. Output in {output_dir}")

    return clean_reads


def run_assembly(
    client: docker.DockerClient,
    sample_id: str,
    clean_reads: dict,
    output_dir: Path,
) -> Path:
    image = DOCKER_IMAGES["assembly"]
    output_dir.mkdir(parents=True, exist_ok=True)

    paired = "r2" in clean_reads
    r1_rel = clean_reads["r1"].relative_to(output_dir).as_posix()

    if paired:
        r2_rel = clean_reads["r2"].relative_to(output_dir).as_posix()
        container_cmd = [
            "--isolate",
            "-1", f"/output/{r1_rel}",
            "-2", f"/output/{r2_rel}",
            "-o", "/output/assembly",
            "--threads", "4",
        ]
    else:
        container_cmd = [
            "--isolate",
            "-s", f"/output/{r1_rel}",
            "-o", "/output/assembly",
            "--threads", "4",
        ]

    logger.info(
        f"Running assembly step (SPAdes --isolate) with image '{image}' for sample '{sample_id}' "
        f"({'paired-end' if paired else 'single-end'})..."
    )

    volumes = {
        to_docker_path(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
    }

    container = client.containers.run(
        image,
        command=container_cmd,
        volumes=volumes,
        detach=True,
    )
    try:
        for line in container.logs(stream=True):
            print(line.decode("utf-8", errors="replace"), end="")
        exit_code = container.wait()["StatusCode"]
    finally:
        container.remove()

    if exit_code != 0:
        logger.error(f"SPAdes failed for sample '{sample_id}' (exit code {exit_code})")
        sys.exit(1)

    contigs = output_dir / "assembly" / "contigs.fasta"
    logger.success(f"Assembly step finished. Contigs written to {contigs}")

    return contigs


def run_quast(
    client: docker.DockerClient,
    sample_id: str,
    contigs: Path,
    output_dir: Path,
) -> None:
    image = DOCKER_IMAGES["quast"]
    output_dir.mkdir(parents=True, exist_ok=True)

    contigs_rel = contigs.relative_to(output_dir).as_posix()

    container_cmd = [
        f"/output/{contigs_rel}",
        "-o", "/output/quast",
        "--threads", "4",
    ]

    logger.info(f"Running QUAST with image '{image}' for sample '{sample_id}'...")

    volumes = {
        to_docker_path(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
    }

    logs = client.containers.run(
        image,
        command=container_cmd,
        volumes=volumes,
        remove=True,
        stdout=True,
        stderr=True,
    )
    logger.info(logs.decode("utf-8", errors="replace"))
    logger.success(f"QUAST finished. Report in {output_dir / 'quast'}")


def run_checkm2(
    client: docker.DockerClient,
    sample_id: str,
    contigs: Path,
    output_dir: Path,
    db_path: Path,
) -> None:
    image = DOCKER_IMAGES["checkm2"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # CheckM2 expects a directory of genome fasta files, not a single file
    # path. For an isolate, that directory just contains the one assembly.
    genome_dir = output_dir / "checkm2_input"
    genome_dir.mkdir(parents=True, exist_ok=True)
    genome_link = genome_dir / "contigs.fasta"
    if not genome_link.exists():
        genome_link.write_bytes(contigs.read_bytes())

    container_cmd = [
        "predict",
        "--input", "/input",
        "--output-directory", "/output/checkm2",
        "--database_path", "/db/CheckM2_database/uniref100.KO.1.dmnd",
        "-x", "fasta",
        "--force",
        "-t", "4",
    ]

    logger.info(f"Running CheckM2 with image '{image}' for sample '{sample_id}'...")

    volumes = {
        to_docker_path(genome_dir.resolve()): {"bind": "/input", "mode": "ro"},
        to_docker_path(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
        to_docker_path(db_path.resolve()): {"bind": "/db", "mode": "ro"},
    }

    container = client.containers.run(
        image,
        command=container_cmd,
        volumes=volumes,
        detach=True,
    )
    try:
        for line in container.logs(stream=True):
            print(line.decode("utf-8", errors="replace"), end="")
        exit_code = container.wait()["StatusCode"]
    finally:
        container.remove()

    if exit_code != 0:
        logger.error(f"CheckM2 failed for sample '{sample_id}' (exit code {exit_code})")
        sys.exit(1)

    logger.success(f"CheckM2 finished. Report in {output_dir / 'checkm2'}")


def run_biosurfdb_search(
    client: docker.DockerClient,
    sample_id: str,
    contigs: Path,
    output_dir: Path,
    db_path: Path,
) -> None:
    """Predict genes on the assembled genome (Prodigal) and search them
    against the local BioSurfDB DIAMOND database.
    """
    image = DOCKER_IMAGES["biosurfdb"]
    biosurfdb_dir = output_dir / "biosurfdb"
    biosurfdb_dir.mkdir(parents=True, exist_ok=True)

    contigs_rel = contigs.relative_to(output_dir).as_posix()

    script = (
        "set -e && "
        f"prodigal -i /output/{contigs_rel} -a /output/biosurfdb/{sample_id}.faa -p single -q && "
        f"diamond blastp "
        f"  -q /output/biosurfdb/{sample_id}.faa "
        f"  -d /db/biosurfdb.dmnd "
        f"  -o /output/biosurfdb/{sample_id}_hits.tsv "
        f"  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle "
        f"  --evalue 1e-5 --max-target-seqs 1 --threads 4"
    )

    logger.info(f"Running BioSurfDB gene prediction + DIAMOND search with image '{image}' for sample '{sample_id}'...")

    volumes = {
        to_docker_path(output_dir.resolve()): {"bind": "/output", "mode": "rw"},
        to_docker_path(db_path.resolve()): {"bind": "/db", "mode": "ro"},
    }

    container = client.containers.run(
        image,
        command=["-c", script],
        volumes=volumes,
        detach=True,
    )
    try:
        for line in container.logs(stream=True):
            print(line.decode("utf-8", errors="replace"), end="")
        exit_code = container.wait()["StatusCode"]
    finally:
        container.remove()

    if exit_code != 0:
        logger.error(f"BioSurfDB search failed for sample '{sample_id}' (exit code {exit_code})")
        sys.exit(1)

    logger.success(f"BioSurfDB search finished. Hits in {biosurfdb_dir / f'{sample_id}_hits.tsv'}")


UNCLASSIFIED = "unclassified / outside surfactant biosynthesis"
BIOSURFDB_FILES = ("biosurfdb.dmnd", "acc2biosurfdb.map", "biosurfdb.map", "biosurfdb.tre")
REPORT_DB_FILES = ("acc2biosurfdb.map", "biosurfdb.map", "biosurfdb.tre")


def missing_biosurfdb_files(db_path: Path, names: tuple[str, ...] = BIOSURFDB_FILES) -> list[str]:
    """Return the names of required BioSurfDB files missing from db_path."""
    return [name for name in names if not (db_path / name).exists()]


def load_two_column_map(path: Path) -> dict[str, str]:
    """Load a tab-separated two-column map file into a dict."""
    mapping: dict[str, str] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping


def parse_newick(text: str) -> dict[str, str]:
    """Parse a simple Newick tree (no branch lengths, numeric labels only)
    into a {child_id: parent_id} map covering every labelled node.
    """
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]

    pos = [0]
    parent: dict[str, str] = {}

    def parse_name() -> str:
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",()":
            pos[0] += 1
        return s[start:pos[0]]

    def parse_clade() -> str:
        if s[pos[0]] == "(":
            pos[0] += 1
            children = []
            while True:
                children.append(parse_clade())
                if s[pos[0]] == ",":
                    pos[0] += 1
                elif s[pos[0]] == ")":
                    pos[0] += 1
                    break
            name = parse_name()
            for child in children:
                parent[child] = name
            return name
        return parse_name()

    parse_clade()
    return parent


def ancestor_chain(node_id: str, parent: dict[str, str], root_id: str) -> list[str] | None:
    """Return the path [node, ..., class] from node_id up to the direct child
    of root_id (the broad class), or None if node_id is not under root_id.
    """
    chain = [node_id]
    current = node_id
    while current in parent:
        p = parent[current]
        if p == root_id:
            return chain
        chain.append(p)
        current = p
    return None


def barh_chart(series, title: str, xlabel: str, path: Path, width: float = 9) -> None:
    """Save a horizontal bar chart from a pandas Series (largest value on top)."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(width, max(4, 0.4 * len(series))))
    ax.barh(
        series.index[::-1],
        series.values[::-1],
        color="#4575b4",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.success(f"Chart written to {path}")


def generate_biosurfdb_report(
    sample_id: str,
    hits_path: Path,
    db_path: Path,
    report_dir: Path,
    root_id: str = "2",
) -> None:
    """Build BioSurfDB reports from a DIAMOND hits TSV (outfmt 6: qseqid sseqid
    pident length evalue bitscore stitle) at three levels of the hierarchy:
    specific categories, subclasses (directly below each broad class) and
    broad classes (direct children of root_id in biosurfdb.tre; default 2,
    "Surfactants"). Each level gets a full CSV table and a bar chart.

    Hits outside root_id are excluded from the summaries, which are
    normalized to the remaining (classified) hits; hits_annotated.csv keeps
    every hit. Pure Python — runs on the host, no container needed.
    """
    import pandas as pd

    if not hits_path.exists() or hits_path.stat().st_size == 0:
        logger.warning(f"No BioSurfDB hits found at {hits_path} — skipping report.")
        return

    report_dir.mkdir(parents=True, exist_ok=True)

    acc2id = load_two_column_map(db_path / "acc2biosurfdb.map")
    id2name = load_two_column_map(db_path / "biosurfdb.map")
    parent_map = parse_newick((db_path / "biosurfdb.tre").read_text(encoding="utf-8", errors="replace"))

    df = pd.read_csv(
        hits_path,
        sep="\t",
        header=None,
        names=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"],
    )
    df["category_id"] = df["sseqid"].map(acc2id).fillna("")
    df["category_name"] = df["category_id"].map(id2name).fillna("unknown")

    def levels_for(category_id: str) -> tuple[str, str]:
        """Return (class_id, subclass_id), or ("", "") if the category is not
        under the root node. A category directly under a class is its own
        subclass.
        """
        chain = ancestor_chain(category_id, parent_map, root_id)
        if chain is None:
            return "", ""
        return chain[-1], chain[-2] if len(chain) >= 2 else chain[-1]

    def name_of(node_id: str) -> str:
        if node_id == "":
            return UNCLASSIFIED
        return id2name.get(node_id, f"unknown (id {node_id})")

    pairs = [levels_for(cid) for cid in df["category_id"]]
    df["class_id"] = [p[0] for p in pairs]
    df["subclass_id"] = [p[1] for p in pairs]
    df["class_name"] = df["class_id"].apply(name_of)
    df["subclass_name"] = df["subclass_id"].apply(name_of)

    # Per-hit annotated table (all hits, all hierarchy levels)
    annotated_path = report_dir / "hits_annotated.csv"
    df.to_csv(annotated_path, index=False)
    logger.success(f"Annotated hits written to {annotated_path}")

    # Summaries exclude hits outside the root node and are normalized to
    # the remaining (classified) hits.
    df_cls = df[df["class_id"] != ""]
    total_hits = len(df_cls)
    logger.info(f"Excluded {len(df) - total_hits} unclassified hits; {total_hits} hits used for the summaries.")
    if total_hits == 0:
        logger.warning("No classified hits left; skipping summaries.")
        return

    # Level 1: specific categories (all)
    cat_counts = df_cls.groupby("category_name").size().sort_values(ascending=False)
    cat_pct = (cat_counts / total_hits * 100).round(2)
    cat_table = cat_counts.reset_index(name="hit_count")
    cat_table["percentage_of_classified_hits"] = cat_pct.values
    cat_table_path = report_dir / "categories_summary.csv"
    cat_table.to_csv(cat_table_path, index=False)
    logger.success(f"Category table written to {cat_table_path}")
    barh_chart(
        cat_pct,
        f"BioSurfDB Functional Categories — Sample {sample_id}",
        "Percentage of classified hits (%)",
        report_dir / "categories_summary.png",
    )

    # Level 2: subclasses (all, with their parent class)
    sub_table = (
        df_cls.groupby(["class_name", "subclass_name"])
        .size()
        .sort_values(ascending=False)
        .reset_index(name="hit_count")
    )
    sub_table["percentage_of_classified_hits"] = (sub_table["hit_count"] / total_hits * 100).round(2)
    sub_table_path = report_dir / "subclasses_summary.csv"
    sub_table.to_csv(sub_table_path, index=False)
    logger.success(f"Subclass table written to {sub_table_path}")
    sub_labels = (sub_table["subclass_name"] + "  (" + sub_table["class_name"] + ")").tolist()
    barh_chart(
        pd.Series(sub_table["percentage_of_classified_hits"].values, index=sub_labels),
        f"BioSurfDB Subclasses — Sample {sample_id}",
        "Percentage of classified hits (%)",
        report_dir / "subclasses_summary.png",
        width=11,
    )

    # Level 3: broad classes (all)
    class_counts = df_cls.groupby("class_name").size().sort_values(ascending=False)
    class_pct = (class_counts / total_hits * 100).round(2)
    class_table = class_counts.reset_index(name="hit_count")
    class_table["percentage_of_classified_hits"] = class_pct.values
    class_table_path = report_dir / "classes_summary.csv"
    class_table.to_csv(class_table_path, index=False)
    logger.success(f"Class-level table written to {class_table_path}")
    barh_chart(
        class_pct,
        f"BioSurfDB Broad Classes — Sample {sample_id}",
        "Percentage of classified hits (%)",
        report_dir / "classes_summary.png",
    )


@click.command()
@click.option(
    "--sample-id",
    required=True,
    help="Sample identifier (defined in config/samples.yaml, unless --hits is used).",
)
@click.option(
    "--config",
    "config_path",
    default=Path(__file__).resolve().parent.parent / "config" / "samples.yaml",
    type=click.Path(path_type=Path),
    help="Path to the samples YAML configuration file.",
)
@click.option(
    "--hits",
    "hits_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help=(
        "Skip the pipeline and only build the BioSurfDB report from an existing "
        "DIAMOND hits TSV (e.g. a fosmid/plasmid assembly searched outside the "
        "pipeline). Output goes to results/<sample-id>/biosurfdb/report."
    ),
)
def main(sample_id: str, config_path: Path, hits_path: Path | None) -> None:
    project_root = Path(__file__).resolve().parent.parent
    biosurfdb_db = project_root / "data" / "biosurfdb"
    checkm2_db = project_root / "data" / "checkm2_db"
    output_dir = project_root / "results" / sample_id

    # Report-only mode: no config lookup, no Docker
    if hits_path is not None:
        missing = missing_biosurfdb_files(biosurfdb_db, REPORT_DB_FILES)
        if missing:
            logger.error(f"Missing BioSurfDB files in {biosurfdb_db}: {', '.join(missing)}")
            sys.exit(1)
        generate_biosurfdb_report(
            sample_id, hits_path, biosurfdb_db, output_dir / "biosurfdb" / "report"
        )
        return

    if not config_path.exists():
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)

    logger.info(f"Loading configuration from: {config_path}")
    samples = load_samples(config_path)

    if sample_id not in samples:
        logger.error(f"Sample '{sample_id}' not found in {config_path}")
        sys.exit(1)

    r1 = project_root / samples[sample_id]["r1"]
    r2 = project_root / samples[sample_id]["r2"] if samples[sample_id].get("r2") else None
    logger.info(f"Sample: {sample_id}")

    # Fail early if a required database is missing, before any long step runs
    if not (checkm2_db / "CheckM2_database" / "uniref100.KO.1.dmnd").exists():
        logger.error(
            f"CheckM2 database not found at {checkm2_db}. "
            f"Download it first with the isolados-biosurf/checkm2 image "
            f"('checkm2 database --download --path /db')."
        )
        sys.exit(1)
    missing = missing_biosurfdb_files(biosurfdb_db)
    if missing:
        logger.error(f"Missing BioSurfDB files in {biosurfdb_db}: {', '.join(missing)}")
        sys.exit(1)

    logger.info("Checking Docker access...")
    client = check_docker()
    if client is None:
        sys.exit(1)
    logger.success("Docker is reachable.")

    clean_reads = run_qc(client, sample_id, r1, r2, output_dir / "qc")
    contigs = run_assembly(client, sample_id, clean_reads, output_dir)
    run_quast(client, sample_id, contigs, output_dir)
    run_checkm2(client, sample_id, contigs, output_dir, checkm2_db)
    run_biosurfdb_search(client, sample_id, contigs, output_dir, biosurfdb_db)
    generate_biosurfdb_report(
        sample_id,
        output_dir / "biosurfdb" / f"{sample_id}_hits.tsv",
        biosurfdb_db,
        output_dir / "biosurfdb" / "report",
    )


if __name__ == "__main__":
    main()
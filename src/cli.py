"""Orchestrator for the isolate genomics pipeline (isolados-biosurf).

Currently implements:
  1. config loading (config/samples.yaml);
  2. Docker daemon check;
  3. QC step (fastp), run as a container. Paired-end only.
  4. assembly step (SPAdes --isolate), run as a container.
  5. QUAST: assembly quality metrics (N50, contig count, etc.).
  6. CheckM2: completeness/contamination estimate for the assembled genome.
  7. BioSurfDB search: Prodigal gene prediction on the assembled genome +
     DIAMOND search against a local BioSurfDB database, summarized by
     biosurfactant pathway category (top-20 table + bar chart).
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
    r2: Path,
    output_dir: Path,
) -> dict:
    image = DOCKER_IMAGES["qc"]
    output_dir.mkdir(parents=True, exist_ok=True)

    clean_r1_name = f"{sample_id}_R1.clean.fastq.gz"
    clean_r2_name = f"{sample_id}_R2.clean.fastq.gz"

    container_cmd = [
        "-i", f"/input/{r1.name}",
        "-I", f"/input/{r2.name}",
        "-o", f"/output/{clean_r1_name}",
        "-O", f"/output/{clean_r2_name}",
        "-j", "/output/fastp.json",
        "-h", "/output/fastp.html",
    ]

    logger.info(f"Running QC step (fastp) with image '{image}' for sample '{sample_id}'...")

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

    return {
        "r1": output_dir / clean_r1_name,
        "r2": output_dir / clean_r2_name,
    }


def run_assembly(
    client: docker.DockerClient,
    sample_id: str,
    clean_reads: dict,
    output_dir: Path,
) -> Path:
    image = DOCKER_IMAGES["assembly"]
    output_dir.mkdir(parents=True, exist_ok=True)

    r1_rel = clean_reads["r1"].relative_to(output_dir).as_posix()
    r2_rel = clean_reads["r2"].relative_to(output_dir).as_posix()

    container_cmd = [
        "--isolate",
        "-1", f"/output/{r1_rel}",
        "-2", f"/output/{r2_rel}",
        "-o", "/output/assembly",
        "--threads", "4",
    ]

    logger.info(f"Running assembly step (SPAdes --isolate) with image '{image}' for sample '{sample_id}'...")

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


def generate_biosurfdb_report(sample_id: str, output_dir: Path, db_path: Path) -> None:
    """Map DIAMOND hits to BioSurfDB functional categories and produce a
    top-20 table plus a bar chart, at the whole-genome level (isolate =
    single genome, no per-bin breakdown needed).
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    hits_path = output_dir / "biosurfdb" / f"{sample_id}_hits.tsv"
    report_dir = output_dir / "biosurfdb" / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    if not hits_path.exists() or hits_path.stat().st_size == 0:
        logger.warning(f"No BioSurfDB hits found at {hits_path} — skipping report.")
        return

    acc2id: dict[str, str] = {}
    with open(db_path / "acc2biosurfdb.map", "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                acc2id[parts[0]] = parts[1]

    id2name: dict[str, str] = {}
    with open(db_path / "biosurfdb.map", "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                id2name[parts[0]] = parts[1]

    df = pd.read_csv(
        hits_path,
        sep="\t",
        header=None,
        names=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"],
    )
    df["category_id"] = df["sseqid"].map(acc2id).fillna("")
    df["category_name"] = df["category_id"].map(id2name).fillna("unknown")

    total_hits = len(df)
    counts = df.groupby("category_name").size().sort_values(ascending=False)
    top20 = counts.head(20)
    top20_pct = (top20 / total_hits * 100).round(2)

    table = top20.reset_index(name="hit_count")
    table["percentage_of_total_hits"] = top20_pct.values
    table_path = report_dir / "top20_categories.csv"
    table.to_csv(table_path, index=False)
    logger.success(f"Top-20 category table written to {table_path}")

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(
        top20_pct.index[::-1],
        top20_pct.values[::-1],
        color="#4575b4",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_xlabel("Percentage of total hits (%)", fontsize=11)
    ax.set_ylabel("")
    ax.set_title(
        f"BioSurfDB Functional Categories — Sample {sample_id}\n"
        "Top 20 categories",
        fontsize=13,
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    chart_path = report_dir / "top20_categories.png"
    fig.savefig(chart_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.success(f"Chart written to {chart_path}")


@click.command()
@click.option(
    "--sample-id",
    required=True,
    help="Sample identifier, as defined in config/samples.yaml.",
)
@click.option(
    "--config",
    "config_path",
    default=Path(__file__).resolve().parent.parent / "config" / "samples.yaml",
    type=click.Path(exists=True, path_type=Path),
    help="Path to the samples YAML configuration file.",
)
def main(sample_id: str, config_path: Path) -> None:
    project_root = Path(__file__).resolve().parent.parent

    logger.info(f"Loading configuration from: {config_path}")
    samples = load_samples(config_path)

    if sample_id not in samples:
        logger.error(f"Sample '{sample_id}' not found in {config_path}")
        sys.exit(1)

    r1 = project_root / samples[sample_id]["r1"]
    r2 = project_root / samples[sample_id]["r2"]
    logger.info(f"Sample: {sample_id}")

    logger.info("Checking Docker access...")
    client = check_docker()
    if client is None:
        sys.exit(1)
    logger.success("Docker is reachable.")

    output_dir = project_root / "results" / sample_id

    clean_reads = run_qc(client, sample_id, r1, r2, output_dir / "qc")
    contigs = run_assembly(client, sample_id, clean_reads, output_dir)
    run_quast(client, sample_id, contigs, output_dir)

    checkm2_db = project_root / "data" / "checkm2_db"
    if not (checkm2_db / "CheckM2_database" / "uniref100.KO.1.dmnd").exists():
        logger.error(
            f"CheckM2 database not found at {checkm2_db}. "
            f"Run 'checkm2 database --download --path {checkm2_db}' first "
            f"(via the isolados-biosurf/checkm2 image)."
        )
        sys.exit(1)
    run_checkm2(client, sample_id, contigs, output_dir, checkm2_db)

    biosurfdb_db = project_root / "data" / "biosurfdb"
    if not (biosurfdb_db / "biosurfdb.dmnd").exists():
        logger.error(
            f"BioSurfDB database not found at {biosurfdb_db}. "
            f"Place biosurfdb.dmnd, acc2biosurfdb.map and biosurfdb.map there first."
        )
        sys.exit(1)
    run_biosurfdb_search(client, sample_id, contigs, output_dir, biosurfdb_db)
    generate_biosurfdb_report(sample_id, output_dir, biosurfdb_db)


if __name__ == "__main__":
    main()
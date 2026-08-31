"""Orchestrator for the isolate genomics pipeline (isolados-biosurf).

Currently implements:
  1. config loading (config/samples.yaml);
  2. Docker daemon check;
  3. QC step (fastp), run as a container. Paired-end only.
  4. assembly step (SPAdes --isolate), run as a container.
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
    run_assembly(client, sample_id, clean_reads, output_dir)


if __name__ == "__main__":
    main()
# IsoBiosurfPipeline 
Automated, containerized pipeline for the genomic characterization of bacterial and archaeal isolates, with a focus on identifying biosurfactant-related functional genes via [BioSurfDB](http://www.biosurfdb.org/).

Sibling project to [metagen-biosurf](https://github.com/ianalvess/metagen-biosurf), adapted from metagenomic community analysis to single-organism isolate genomes.

## Table of contents

- [Overview](#overview)
- [Pipeline](#pipeline)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)
- [Project structure](#project-structure)
- [Scientific background](#scientific-background)

## Overview

Biosurfactant screening depends on functional gene content rather than taxonomy alone, since biosynthesis pathways can be encoded across diverse and sometimes unclassified taxa. This pipeline processes Illumina paired-end reads from a single isolate through quality control, assembly, quality assessment, and functional annotation, producing a ranked summary of biosurfactant-related gene categories detected in the genome.

Each pipeline stage runs in its own Docker container, orchestrated by a Python CLI. No cluster scheduler (Snakemake, Nextflow) is used by design — the orchestration logic lives entirely in `src/cli.py`.

## Pipeline

| Stage | Tool | Purpose |
|---|---|---|
| QC | fastp | Adapter/quality trimming of raw reads |
| Assembly | SPAdes (`--isolate`) | De novo genome assembly |
| Quality assessment | QUAST | Assembly contiguity metrics (N50, contig count, genome size) |
| Quality assessment | CheckM2 | Genome completeness and contamination estimate |
| Functional annotation | Prodigal + DIAMOND | Gene prediction and search against BioSurfDB |

## Requirements

- Docker Desktop (with the WSL2 backend, on Windows)
- Python 3.10+
- ~10GB free disk space for reference databases

## Installation

**1. Clone and build the Docker images:**

```bash
git clone https://github.com/ianalvess/isolados-biosurf.git
cd isolados-biosurf

docker build -t isolados-biosurf/fastp docker/fastp
docker build -t isolados-biosurf/spades docker/spades
docker build -t isolados-biosurf/quast docker/quast
docker build -t isolados-biosurf/checkm2 docker/checkm2
docker build -t isolados-biosurf/biosurfdb docker/biosurfdb
```

**2. Download reference databases:**

CheckM2 (~3.5GB):

```bash
docker run --rm -v <absolute-path-to>/data/checkm2_db:/db isolados-biosurf/checkm2 database --download --path /db
```

BioSurfDB: place `biosurfdb.dmnd`, `acc2biosurfdb.map`, and `biosurfdb.map` in `data/biosurfdb/`.

**3. Install Python dependencies:**

```bash
pip install click docker loguru pyyaml pandas matplotlib
```

## Configuration

Register each sample's raw read paths in `config/samples.yaml`:

```yaml
samples:
  M11:
    r1: "data/M11/M11_R1.fq.gz"
    r2: "data/M11/M11_R2.fq.gz"
```

Reads must be placed inside the project's `data/` directory.

## Usage

```bash
python src/cli.py --sample-id M11
```

Runs the full pipeline end to end for the specified sample.

## Output

Results are written to `results/<sample_id>/`:

```
results/<sample_id>/
├── qc/                     # clean reads + fastp HTML/JSON report
├── assembly/                # contigs.fasta + SPAdes logs
├── quast/                   # assembly quality report
├── checkm2/                 # completeness/contamination report
└── biosurfdb/
    ├── <sample_id>_hits.tsv # raw DIAMOND hits
    └── report/
        ├── top20_categories.csv
        └── top20_categories.png
```

## Project structure

```
isolados-biosurf/
├── config/
│   └── samples.yaml
├── data/
│   ├── <sample_id>/         # raw reads
│   ├── checkm2_db/
│   └── biosurfdb/
├── docker/
│   ├── fastp/
│   ├── spades/
│   ├── quast/
│   ├── checkm2/
│   └── biosurfdb/
├── results/
├── src/
│   └── cli.py
└── logs/
```

## Scientific background

Biosurfactants are microbially produced surface-active compounds with applications in bioremediation, enhanced oil recovery, and green chemistry. This pipeline targets pure, single-organism isolates rather than environmental metagenomes: working with an isolate genome removes the binning and reconciliation uncertainty inherent to metagenomic assembly, giving an unambiguous gene inventory for the organism in question.

Assembly quality (QUAST) and genome completeness/contamination (CheckM2) are assessed before functional annotation, ensuring BioSurfDB hits are attributed to a reliable genome rather than assembly artifacts. The resulting candidate biosurfactant-related genes can then guide downstream work such as pathway confirmation, culture-based validation, or comparison against isolates from other extreme environments.
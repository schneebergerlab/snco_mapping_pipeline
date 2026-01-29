<picture>
  <source
    srcset="images/logo.png"
    media="(orientation: portrait)" />
  <img src="images/logo.png" alt="" width=600 />
</picture>

A snakemake pipeline for mapping reads from various types of plate or droplet-based single-cell platforms, for crossover/haplotyping analysis using `coelsch`. The main coelsch tool repo is available [here](https://www.github.com/schneebergerlab/coelsch.git).

Several droplet-based single-cell modalities are currently supported, including 10x 3' RNA v3 and v4, as well as BD Rhapsody 3' RNA and 10x ATAC. For these datasets one or more pairs of fastq files containing reads from all barcodes should be supplied.

Plate-based methods like Takara PicoPlex WGA or analysis of whole-genome resequencing datasets for backcross/F2 populations (see [Rowan et al. 2015](https://doi.org/10.1534/g3.114.016501) for an example) are supported using individual input fastq files for each barcode and/or individual.

## Installation

    git clone https://github.com/schneebergerlab/coelsch_mapping_pipeline.git
    cd coelsch_mapping_pipeline
    pip install -e .

## Usage

### Initialise config

    coelsch_pipeline init

This creates a `coelsch_mapping_config.yaml` in your current directory, prompting you to fill in some of the required fields. The pipeline run requires a directory structure specified by the config file. The `annotation_dir` directory should contain the reference genomes and annotations, plus optional predefined variants, in fasta, gtf, and vcf formats respectively. The `raw_data_dir` directory should contain the input fastq files. The `results_dir` directory will be populated by the aligned data and haplotyping results during the pipeline run.

### Run the pipeline

    coelsch_pipeline run

The snakefiles directing the pipeline structure are located with the installation of the tool and are not present in the output directory. This is to make it easier to create multiple independent runs of the pipeline without copying/symlinking snakefiles.

The `coelsch_pipeline` cli is a fairly straightforward wrapper of the snakemake cli. You can supply any Snakemake arguments (e.g. `--profile`, `--use-conda`, etc.) which will be passed directly onwards to snakemake:

    coelsch_pipeline run --cores 4 --use-conda --rerun-incomplete

The only arguments that should be avoided are `--snakefile` and `--configfile` since these are managed by the cli.

### Requirements & conda support

The pipeline itself has few requirements - just `snakemake`, `jinja2`, and `click`. None of the tools used for the analysis are installed by the `setup.py` script itself, but are instead managed by snakemake using conda. The default environment yamls used by snakemake are stored in the subdirectory `coelsch_mapping_pipeline/coelsch_pipeline/pipeline/env_yamls/` of this repo. If you have problems with the conda installations of these environments, alternative yamls or existing conda environment names can be provided in the `conda_envs` section of the `coelsch_mapping_config.yaml` file.

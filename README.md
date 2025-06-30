# snco_pipeline

Pipeline tool for mapping reads for crossover analysis with snco

## Installation

    pip install -e .

## Usage

### Initialise config

    snco_pipeline init-config

This creates a `config.yaml` in your current directory, prompting you to fill in required fields.

### Run the pipeline

    snco_pipeline run

You can pass any Snakemake arguments (e.g. `--profile`, `--use-conda`, etc.):

    snco_pipeline run --cores 4 --use-conda --rerun-incomplete

You can also override config values inline:

    snco_pipeline run --set outfile=test.txt

## Example output

The default pipeline just touches an output file defined in the config.
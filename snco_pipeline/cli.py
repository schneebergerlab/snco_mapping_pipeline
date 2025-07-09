import os
import shutil
import sys
from pathlib import Path
import yaml

import click
from jinja2 import Environment, FileSystemLoader, meta


BASE_DIR = Path(__file__).parent
PIPELINE_DIR = BASE_DIR / "pipeline"
CONFIG_DIR = PIPELINE_DIR / "config"
ANNOTATIONS_DIR = PIPELINE_DIR / "annotations"
SNAKEFILE = PIPELINE_DIR / "Snakefile"

CONFIG_VARIABLES = [
    'annotation_dir', 'raw_data_dir', 'results_dir',
    'single_cell_method', 'reference_genotype',
    'dataset_name', 'dataset_genotype_name', 'parent1_name', 'parent2_name'
]

BARCODE_WHITELIST_FNS = {
    '10x_rna_v4': ['3M-3pgex-may-2023.txt',],
    '10x_rna_v3': ['3M-february-2018.txt',],
    '10x_atac': ['737K-cratac-v1.txt',],
    'bd_rna': ['BD_CLS1.txt', 'BD_CLS2.txt', 'BD_CLS3.txt'],
    'takara_dna': None,
    'plate_wgs': None,
}

if not SNAKEFILE.exists():
    click.echo(f"Error: Snakefile not found at {SNAKEFILE}", err=True)
    sys.exit(1)


@click.group()
def cli():
    """CLI for the pipeline."""
    pass


@cli.command("init")
@click.argument("destination", type=click.Path(), default='snco_mapping_config.yaml')
@click.option("-f", "--force", is_flag=True, default=False)
def init_config(destination, force):
    """Prompt for missing values and render config.yaml from template."""
    dest = Path(destination)
    if dest.exists() and not force:
        click.secho(f"{dest} already exists. Refusing to overwrite.", err=True, fg="red")
        sys.exit(1)

    env = Environment(loader=FileSystemLoader(CONFIG_DIR))
    template_name = "default_config.yaml.j2"
    template_source = env.loader.get_source(env, template_name)[0]

    # Find undeclared variables in the template
    parsed = env.parse(template_source)
    variables = meta.find_undeclared_variables(parsed)

    # Prompt for each variable
    context = {}
    for var in CONFIG_VARIABLES:
        value = click.prompt(click.style(f"Enter value for '{var}'", fg='bright_yellow'), default='')
        if value != '':
            context[var] = value

    tech_type = context.get('single_cell_method', '10x_rna_v4')
    context['barcode_whitelist_fns'] = str(BARCODE_WHITELIST_FNS.get(tech_type, None))

    # Render and write output
    template = env.get_template(template_name)
    rendered = template.render(context)
    dest.write_text(rendered)
    click.secho(f"Config written to {dest}", fg="green")
    config = yaml.safe_load(rendered)
    for directory in ['annotation_dir', 'raw_data_dir', 'results_dir']:
        dir_name = config[directory]
        if not os.path.exists(dir_name):
            click.secho(f'Creating {directory}: {dir_name}', fg="bright_yellow")
            os.mkdir(dir_name)


@cli.command("run", context_settings={"ignore_unknown_options": True, "allow_extra_args": True})
@click.option("--configfile", "-c", type=click.Path(exists=True), default="snco_mapping_config.yaml",
              help="Path to config file")
@click.pass_context
def run_pipeline(ctx, configfile):
    """Run the pipeline."""

    args = [
        "--snakefile", str(SNAKEFILE.resolve()),
        "--configfile", str(configfile),
    ]

    # Add any unknown extra CLI args
    args += ctx.args

    click.secho(f"Running: snakemake {' '.join(args)}", fg="bright_yellow")

    try:
        # Snakemake ≥7.x
        from snakemake.cli import main as snakemake
    except ImportError:
        # Snakemake ≤6.x
        from snakemake import main as snakemake

    snakemake(args)

if __name__ == "__main__":
    cli()

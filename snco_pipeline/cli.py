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
}

if not SNAKEFILE.exists():
    click.echo(f"Error: Snakefile not found at {SNAKEFILE}", err=True)
    sys.exit(1)


@click.group()
def cli():
    """CLI for the pipeline."""
    pass


@cli.command("init")
@click.argument("destination", type=click.Path(), default='config.yaml')
@click.option("-f", "--force", is_flag=True, default=False)
def init_config(destination, force):
    """Prompt for missing values and render config.yaml from template."""
    dest = Path(destination)
    if dest.exists() and not force:
        click.echo(f"{dest} already exists. Refusing to overwrite.", err=True)
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
        value = click.prompt(f"Enter value for '{var}'", default='')
        if value != '':
            context[var] = value

    tech_type = context.get('single_cell_method', '10x_rna_v4')
    context['barcode_whitelist_fns'] = str(BARCODE_WHITELIST_FNS[tech_type])

    # Render and write output
    template = env.get_template(template_name)
    rendered = template.render(context)
    dest.write_text(rendered)
    click.echo(f"Config written to {dest}")
    config = yaml.safe_load(rendered)
    for directory in ['annotation_dir', 'raw_data_dir', 'results_dir']:
        dir_name = config[directory]
        if not os.path.exists(dir_name):
            click.echo(f'Creating {directory}: {dir_name}')
            os.mkdir(dir_name)
    annot_dir = config['annotation_dir']
    barcode_whitelist_fns = [ANNOTATIONS_DIR / fn for fn in BARCODE_WHITELIST_FNS[tech_type]]
    for fn in barcode_whitelist_fns:
        click.echo(f'copying barcode file {fn} to {annot_dir}')
        shutil.copy(fn, annot_dir)



@cli.command("run", context_settings={"ignore_unknown_options": True, "allow_extra_args": True})
@click.option("--configfile", "-c", type=click.Path(exists=True), default="config.yaml",
              help="Path to config file")
@click.option("--set", "overrides", multiple=True, help="Override config values (key=value)")
@click.pass_context
def run_pipeline(ctx, configfile, overrides):
    """Run the pipeline."""
    # Convert key=value overrides into a flat list for --config
    config_overrides = []
    for item in overrides:
        if "=" not in item:
            click.echo(f"Ignoring malformed --set '{item}'", err=True)
            continue
        config_overrides.append(item)

    args = [
        "--snakefile", str(SNAKEFILE.resolve()),
        "--configfile", str(configfile),
    ]

    if config_overrides:
        args += ["--config"] + config_overrides

    # Add any unknown extra CLI args
    args += ctx.args

    click.echo(f"Running: snakemake {' '.join(args)}")

    try:
        # Snakemake ≥7.x
        from snakemake.cli import main as snakemake
    except ImportError:
        # Snakemake ≤6.x
        from snakemake import main as snakemake

    snakemake(args)

if __name__ == "__main__":
    cli()

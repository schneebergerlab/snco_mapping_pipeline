from setuptools import setup, find_packages

setup(
    name="snco_pipeline",
    version="0.1.1",
    description="Pipeline tool for mapping reads for crossover analysis with snco",
    author="Matthew Parker",
    author_email="mparker@mpipz.mpg.de",
    packages=["snco_pipeline"],
    install_requires=[
        "click>=8.0",
        "jinja2",
        "snakemake",
    ],
    entry_points={
        "console_scripts": [
            "snco_pipeline = snco_pipeline.cli:cli"
        ]
    },
    include_package_data=True,
    zip_safe=False,
)

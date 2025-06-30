import os
from snakemake.io import ancient, expand


def resolve_path(basedir, path):
    if isinstance(path, list):
        return [os.path.join(basedir,  p) for p in path]
    return os.path.join(basedir,  path)


def make_path_getter(config_key):
    def factory(config):
        base = config
        for part in config_key.split('/'):
            base = base[part]
        def _wrapped(path):
            return resolve_path(base, path)
        return _wrapped
    return factory


annotations_getter = make_path_getter("annotation_dir")
raw_data_getter = make_path_getter("raw_data_dir")
results_getter = make_path_getter("results_dir")


def make_annotation_getter(config_key):
    def factory(config):
        base = config
        for part in config_key.split('/'):
            base = base[part]
        annotations = annotations_getter(config)
        def _wrapped(annot_key):
            return ancient(annotations(base[annot_key]))
        return _wrapped
    return factory


fasta_getter = make_annotation_getter("annotations/fasta_fns")
gtf_getter = make_annotation_getter("annotations/gtf_fns")
barcode_whitelist_getter = make_annotation_getter("annotations/barcode_whitelist_fns")

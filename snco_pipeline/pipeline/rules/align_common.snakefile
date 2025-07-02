import os
import math
import gzip
from glob import glob

def get_star_index_input(wc):
    '''input for STAR genome index generation'''
    input_ = {
        'fasta': get_fasta(wc.ref),
        'gtf': get_gtf(wc.ref),
    }

    if wc.ref != wc.qry:
        input_['vcf'] = annotation('vcf/star_consensus/msyd/{geno_group}.{qry}.vcf')
    return input_


def get_read_length(fastq_fn):
    '''get assumed length of reads in a fastq file from the first read'''
    if fastq_fn.endswith('.gz'):
        open_method = gzip.open
    else:
        open_method = open
    with open_method(fastq_fn) as f:
        next(f) # skip first read id
        first_read_seq = next(f).decode().strip()
    return len(first_read_seq)


def get_sj_overhang_size(wc):
    '''get the splice junction overhang size for a genome index from the data'''
    if config['alignment']['star'].get('sjdb_overhang_length') is not None:
        return config['alignment']['star']['sjdb_overhang_length']
    fastq_fns = glob(
        raw_data(f'*{config["file_suffixes"]["read2"]}')
    )[:5] # just check a few random files or this gets slow
    max_overhang = 0
    for fq_fn in fastq_fns:
        max_overhang = max(max_overhang, get_read_length(fq_fn))
    return max_overhang - 1


def calculate_genome_sa_index_nbases(wc, input):
    '''estimate the suffix array nbases for star from the size of the reference genome'''
    fai_path = str(input.fasta) + '.fai'
    if not os.path.exists(fai_path):
        if config['alignment']['star'].get('genome_sa_index_nbases') is not None:
            return config['alignment']['star']['genome_sa_index_nbases']
        try:
            import pysam
            pysam.index(input.fasta)
        except:
            raise ValueError(
                f'Cannot find/create fai file for "{input.fasta}". Please provide either .fai or '
                'specify alignment:star:genome_sa_index_nbases in the config file'
            )
    total_length = 0
    with open(fai_path) as f:
        for line in f:
            _, ln, *_ = line.strip().split('\t')
            total_length += int(ln)
    return min(14, math.floor(math.log2(total_length) / 2 - 1))


rule build_STAR_index:
    """
    Builds a STAR genome index for alignment. Uses either reference genome alone, or in conjunction with a 
    VCF file to produce a "STAR consensus" variant transformed genome.
    """
    input:
        unpack(get_star_index_input)
    output:
        directory(annotation('star_indexes/{geno_group}/{ref}.{qry}.star_idx'))
    threads: 12
    resources:
        mem_mb=12 * 1024,
    params:
        overhang=get_sj_overhang_size,
        genome_sa_index_nbases=calculate_genome_sa_index_nbases,
        vcf_flag=lambda wc, input: f'--genomeTransformVCF {input.vcf}' if hasattr(input, 'vcf') else '',
        transform_flag=lambda wc, input: '--genomeTransformType Haploid' if hasattr(input, 'vcf') else '',
    conda:
        '../env_yamls/star.yaml'
    shell:
        '''
        mkdir {output}
        STAR \
          --outTmpDir _STARtmp_{wildcards.geno_group}.{wildcards.ref}.{wildcards.qry} \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles {input.fasta} \
          --sjdbGTFfile {input.gtf} \
          {params.vcf_flag} \
          {params.transform_flag} \
          --genomeSAindexNbases {params.genome_sa_index_nbases} \
          --sjdbOverhang {params.overhang}
        '''


def get_geno_group(cond):
    '''creates a geno_group wildcard from the config by joining together all the genotypes'''
    dataset = config['datasets'][cond]
    ref = dataset['reference_genotype']
    qry_names = set()
    for geno in dataset['genotypes'].values():
        for qry in geno.values():
            if qry != ref:
                qry_names.add(qry)
    qry_names = '_'.join(sorted(qry_names))
    return f'{ref}_{qry_names}'


def get_star_fastq_input(cond, fastq_type):
    '''use globbing to identify all the input fastqs for a dataset from input_file_basenames'''
    fastq_fn_globs = expand(
        raw_data('{sample_name}{file_suffix}'),
        sample_name=config['datasets'][cond]['input_file_basenames'],
        file_suffix=config['file_suffixes'][fastq_type]
    )
    fastq_fns = []
    for fq_fn_glb in fastq_fn_globs:
        fq_fns = glob(fq_fn_glb)
        if not fq_fns:
            raise ValueError(f'could not identify any files matching pattern "{fq_fn_glb}"')
        fastq_fns += fq_fns
    return fastq_fns

import math
import gzip
from glob import glob

def get_star_index_input(wc):
    input_ = {
        'fasta': get_fasta(wc.ref),
        'gtf': get_gtf(wc.ref),
    }

    if wc.ref != wc.qry:
        input_['vcf'] = annotation('vcf/star_consensus/msyd/{geno_group}.{qry}.vcf')
    return input_


def get_read_length(fastq_fn):
    with gzip.open(fastq_fn) as f:
        next(f) # skip first read id
        first_read_seq = next(f).decode().strip()
    return len(first_read_seq)


def get_sj_overhang_size(wc):
    fastq_fns = glob(
        raw_data(f'*{config["file_suffixes"]["read2"]}')
    )[:5] # just check a few random files or this gets slow
    max_overhang = 0
    for fq_fn in fastq_fns:
        max_overhang = max(max_overhang, get_read_length(fq_fn))
    return max_overhang - 1


def calculate_genome_sa_index_nbases(wc, input):
    fai_path = str(input.fasta) + '.fai'
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

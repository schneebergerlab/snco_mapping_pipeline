from glob import glob

include: './align_common.snakefile'


DATASET_SAMPLE_NAME_MAPPING = {}
for dataset_name in config['datasets']:
    # config only stores dataset_name -> sample_name_glob relationship
    # use globbing to determine the sample_name -> dataset_name relationship
    sample_names = glob_wildcards(
        raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}'),
        files=get_star_fastq_input(dataset_name, 'read1')
    ).sample_name
    for sn in sample_names:
        DATASET_SAMPLE_NAME_MAPPING[sn] = dataset_name


def STAR_consensus_input(wc):
    '''
    full input to STAR consensus (plate mode WITHOUT STARsolo).
    input is:
      - index: the STAR index generated from the reference genome plus VCF transform
      - read, mate: the fastq files for the sample - represents a single individual/barcode
    '''
    dataset_name = DATASET_SAMPLE_NAME_MAPPING[wc.sample_name]
    dataset = config['datasets'][dataset_name]
    tech_type = dataset['technology']
    ref_name = dataset['reference_genotype']
    geno_group = get_geno_group(dataset_name)
    return {
        'index': annotation(f'star_indexes/{geno_group}/{ref_name}.{{qry}}.star_idx'),
        'read': raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}'),
        'mate': raw_data(f'{{sample_name}}{config["file_suffixes"]["read2"]}'),
    }


def get_transform_flag(wc):
    '''STAR genome transform flag - necessary if index has a VCF transformation'''
    dataset_name = DATASET_SAMPLE_NAME_MAPPING[wc.sample_name]
    dataset = config['datasets'][dataset_name]
    ref = dataset['reference_genotype']
    if wc.qry != ref:
        return '--genomeTransformOutput SAM'
    else:
        return ''


def get_temp_dir(wc, output):
    '''
    create a private directory for each dataset_name/query combination, since STAR does not work well 
    when multiple processes are run in the same directory
    '''
    output_dir = os.path.split(output.bam)[0]
    tmp_dir = os.path.join(output_dir, f'{wc.sample_name}.{wc.qry}.tmpdir')
    return tmp_dir


rule STAR_consensus:
    '''
    map reads for a single individual/barcode using STARsolo and a VCF-transformed genome index.
    '''
    input:
        unpack(STAR_consensus_input)
    output:
        bam=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam')),
    params:
        star_tmp_dir=get_temp_dir,
        transform_flag=get_transform_flag,
        filter_multimap_nmax=config['alignment']['star']['filter_multimap_nmax'],
        filter_mismatch_nmax=config['alignment']['star']['filter_mismatch_nmax'],
        align_mates_gap_max=config['alignment']['star']['atac']['mates_gap_max'],
    log:
        progress=results('logs/{sample_name}.{qry}.STAR_progress.log'),
        final=results('logs/{sample_name}.{qry}.STAR_final.log'),
        main=results('logs/{sample_name}.{qry}.STAR.log')
    threads: 6 # usually small files, fewer threads needed than for droplets
    resources:
        mem_mb=lambda wildcards, threads: (threads + 4) * 2048,
    conda:
        get_conda_env('star')
    shell:
        '''
        mkdir -p {params.star_tmp_dir}
        RELPATH=$(realpath --relative-to="{params.star_tmp_dir}" ".")
        cd {params.star_tmp_dir}
        STAR \
          --runThreadN {threads} \
          --genomeDir "${{RELPATH}}/{input.index}" \
          --readFilesCommand "zcat" \
          --readFilesIn "${{RELPATH}}/{input.read}" "${{RELPATH}}/{input.mate}" \
          --alignIntronMax 1 \
          --alignMatesGapMax {params.align_mates_gap_max} \
          --outFilterMultimapNmax {params.filter_multimap_nmax} \
          --outFilterMismatchNmax {params.filter_mismatch_nmax} \
          --outSAMtype "BAM" "Unsorted" \
          --outSAMattrRGline "ID:{wildcards.qry}" \
          --outSAMattributes "NH" "HI" "AS" "nM" "RG" \
          {params.transform_flag}

        cd $RELPATH
        mv {params.star_tmp_dir}/Log.progress.out {log.progress}
        mv {params.star_tmp_dir}/Log.final.out {log.final}
        mv {params.star_tmp_dir}/Log.out {log.main}

        samtools sort \
          -m 1G -@ {threads} \
          -o {output.bam} \
          {params.star_tmp_dir}/Aligned.out.bam
        samtools index {output.bam}
          
        rm -rf {params.star_tmp_dir}
        '''


rule sort_bam_by_name:
    '''
    name-sort the output of star consensus to ensure consistent order of read-ids across haplotypes
    '''
    input:
        bam=results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam'),
        bai=results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam.bai')
    output:
        bam=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'))
    threads: 4
    resources:
        mem_mb=20_000,
    conda:
        get_conda_env('htslib')
    shell:
        '''
        samtools sort -n -@ {threads} {input.bam} > {output.bam}
        '''


def get_merge_haps_input(wc):
    '''Expand all haplotypes that have been aligned to for a sample/barcode to create merge input'''
    dataset_name = DATASET_SAMPLE_NAME_MAPPING[wc.sample_name]
    dataset = config['datasets'][dataset_name]
    qrys = set()
    for geno in dataset['genotypes'].values():
        for qry in geno.values():
            qrys.add(qry)
    return expand(
        results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'),
        sample_name=wc.sample_name, qry=qrys
    )


rule merge_name_sorted_hap_bams:
    '''
    merge all haplotype-specific alignments for a sample/barcode into a single name-sorted bam file for haplotype collapsing
    '''
    input:
        bams=get_merge_haps_input
    output:
        bam=temp(results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')),
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1024,
    conda:
        get_conda_env('htslib')
    shell:
        '''
        samtools merge -@ {threads} -n {output.bam} {input.bams}
        '''


rule collapse_alignments:
    '''
    Uses snco script collapse_ha_specific_alns.py to select the best haplotype alignment(s) for each read
    and outputs a single alignment with a new ha tag that indicates which haplotype(s) is/are best.
    Also adds new RG tag to bam file to indicate the sample name, so that this can be identified later
    after sample-wise merging
    '''
    input:
        bam=results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')
    output:
        bam=results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        bai=results('aligned_data/single_barcodes/{sample_name}.sorted.bam.bai')
    resources:
        mem_mb=5_000,
    threads: 4
    conda:
        get_conda_env('snco')
    shell:
        '''
        collapse_ha_specific_alns.py \
          -o {output.bam}.unsorted.bam {input.bam}
        samtools addreplacerg  \
          -r "ID:{wildcards.sample_name}" \
          {output.bam}.unsorted.bam | \
        samtools sort -@ {threads} \
          -T ${{TMPDIR}}/{wildcards.sample_name} \
          -o {output.bam} -
        # add sample_name as read group to 
        samtools index {output.bam}
        rm {output.bam}.unsorted.bam
        '''


def merge_samples_input(wc):
    '''input for sample/barcode-wise merging'''
    sample_names = glob_wildcards(
        raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}'),
        files=get_star_fastq_input(wc.dataset_name, 'read1')
    ).sample_name
    return expand(
        results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        sample_name=sample_names
    )


rule list_bam:
    '''create a list of bam files representing all the samples for a dataset, for merging'''
    input:
        bams=merge_samples_input
    output:
        txt=temp(results("{dataset_name}.list"))
    run:
        with open(output.txt, 'w') as f:
            f.write('\n'.join(input.bams))


rule merge_bams:
    '''sample-wise barcode merging'''
    input:
        results("{dataset_name}.list")
    output:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai'),
    conda:
        get_conda_env('htslib')
    threads: 24
    shell:
        '''
        ulimit -n 5000
        samtools merge -@ {threads} -b {input} -o {output.bam}
        samtools index {output.bam}
        '''

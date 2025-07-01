from glob import glob

include: './align_common.snakefile'


COND_SAMPLE_NAME_MAPPING = {}
for cond in config['datasets']:
    # config only stores cond -> sample_name_glob relationship
    # use globbing to determine the sample_name -> cond relationship
    sample_names = glob_wildcards(
        raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}'),
        files=get_star_fastq_input(cond, 'read1')
    ).sample_name
    for sn in sample_names:
        COND_SAMPLE_NAME_MAPPING[sn] = cond


def STAR_consensus_input(wc):
    cond = COND_SAMPLE_NAME_MAPPING[wc.sample_name]
    dataset = config['datasets'][cond]
    tech_type = dataset['technology']
    ref_name = dataset['reference_genotype']
    geno_group = get_geno_group(cond)
    return {
        'index': annotation(f'star_indexes/{geno_group}/{ref_name}.{{qry}}.star_idx'),
        'read': raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}'),
        'mate': raw_data(f'{{sample_name}}{config["file_suffixes"]["read2"]}'),
    }


def get_transform_flag(wc):
    cond = COND_SAMPLE_NAME_MAPPING[wc.sample_name]
    dataset = config['datasets'][cond]
    ref = dataset['reference_genotype']
    if wc.qry != ref:
        return '--genomeTransformOutput SAM'
    else:
        return ''


def get_temp_dir(wc, output):
    output_dir = os.path.split(output.bam)[0]
    tmp_dir = os.path.join(output_dir, f'{wc.sample_name}.{wc.qry}.tmpdir')
    return tmp_dir


rule STAR_consensus:
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
        '../env_yamls/star.yaml'
    shell:
        '''
        mkdir -p {params.star_tmp_dir}
        RELPATH=$(realpath --relative-to="{params.star_tmp_dir}" ".")
        cd {params.star_tmp_dir}
        STAR \
          --runThreadN {threads} \
          --genomeDir "${{RELPATH}}/{input.index}" \
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
    input:
        bam=results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.sorted.bam')
    output:
        bam=temp(results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'))
    threads: 4
    resources:
        mem_mb=20_000,
    conda:
        '../env_yamls/snco.yaml'
    shell:
        '''
        samtools sort -n -@ {threads} {input.bam} > {output.bam}
        '''


def get_merge_haps_input(wc):
    cond = COND_SAMPLE_NAME_MAPPING[wc.sample_name]
    dataset = config['datasets'][cond]
    qrys = set()
    for geno in dataset['genotypes'].values():
        for qry in geno.values():
            qrys.add(qry)
    return expand(
        results('aligned_data/single_barcodes/haploid/{sample_name}.{qry}.namesorted.bam'),
        sample_name=wc.sample_name, qry=qrys
    )


rule merge_name_sorted_hap_bams:
    input:
        bams=get_merge_haps_input
    output:
        bam=temp(results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')),
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1024,
    conda:
        '../env_yamls/snco.yaml'
    shell:
        '''
        samtools merge -@ {threads} -n {output.bam} {input.bams}
        '''


rule collapse_alignments:
    input:
        bam=results('aligned_data/single_barcodes/{sample_name}.namesorted.bam')
    output:
        bam=results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        bai=results('aligned_data/single_barcodes/{sample_name}.sorted.bam.bai')
    resources:
        mem_mb=5_000,
    threads: 4
    conda:
        '../env_yamls/snco.yaml'
    shell:
        '''
        collapse_ha_specific_alns.py \
          -o {output}.unsorted.bam {input.bam}
        samtools addreplacerg  \
          -r "ID:{wildcards.sample_name}" \
          {output}.unsorted.bam | \
        samtools sort -@ {threads} \
          -T ${{TMPDIR}}/{wildcards.sample_name} \
          -o {output.bam} -
        # add sample_name as read group to 
        samtools index {output.bam}
        rm {output.bam}.unsorted.bam
        '''


def merge_samples_input(wc):
    sample_names = glob_wildcards(
        raw_data(f'{{sample_name}}{config["file_suffixes"]["read1"]}'),
        files=get_star_fastq_input(wc.cond, 'read1')
    ).sample_name
    return expand(
        results('aligned_data/single_barcodes/{sample_name}.sorted.bam'),
        sample_name=sample_names
    )


rule merge_samples:
    input:
        bams=merge_samples_input
    output:
        bam=results('aligned_data/{cond}.sorted.bam'),
        bai=results('aligned_data/{cond}.sorted.bam.bai')
    resources:
        mem_mb=20_000,
    threads: 12
    conda:
        '../env_yamls/snco.yaml'
    shell:
        '''
        samtools merge -@ {threads} {output.bam} {input.bams}
        samtools index {output.bam}
        '''

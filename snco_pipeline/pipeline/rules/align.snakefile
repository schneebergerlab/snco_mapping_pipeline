import math
import gzip
from glob import glob

def get_star_index_input(wc):
    input_ = {
        'fasta': get_fasta(wc.ref),
        'gtf': get_gtf(wc.ref),
    }

    if wc.ref != wc.qry:
        if len(wc.geno_group.split('_')) > 2:
            input_['vcf'] = annotation('vcf/star_consensus/msyd/{geno_group}.{qry}.vcf')
        else:
            input_['vcf'] = annotation('vcf/star_consensus/syri/{ref}.{qry}.vcf')
    return input_


def get_read_length(fastq_fn):
    with gzip.open(fastq_fn) as f:
        next(f) # skip first read id
        first_read_seq = next(f).decode().strip()
    return len(first_read_seq)


def get_sj_overhang_size(wc):
    fastq_fns = glob(
        raw_data('{sample_name}.2.fastq.gz')
    )
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
    '''Create the index required for alignment with STAR'''
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


rule cell_barcode_rc:
    input:
        barcode=raw_data(f'{{sample_name}}{config["file_suffixes"]["atac"]["barcode"]}'),
    output:
        barcode=raw_data('{sample_name}.bc_rc.fastq.gz'),
    conda:
        '../env_yamls/seqtk.yaml'
    shell:
        '''
        seqtk seq -r {input.barcode} | bgzip > {output.barcode}
        '''


def get_geno_group(wc):
    dataset = config['datasets'][wc.cond]
    ref = dataset['reference_genotype']
    qry_names = set()
    for geno in dataset['genotypes'].values():
        for qry in geno.values():
            if qry != ref:
                qry_names.add(qry)
    qry_names = '_'.join(sorted(qry_names))
    return f'{ref}_{qry_names}'


def get_star_fastq_input(cond, dtype, fastq_type):
    return expand(
        raw_data('{sample_name}{file_suffix}'),
        sample_name=config['datasets'][cond]['input_file_basenames'],
        file_suffix=config['file_suffixes'][dtype][fastq_type]
    )


def STAR_consensus_input(wc):
    dataset = config['datasets'][wc.cond]
    tech_type = dataset['technology']
    ref_name = dataset['reference_genotype']
    geno_group = get_geno_group(wc)
    input_ = {
        'index': annotation(f'star_indexes/{geno_group}/{ref_name}.{{qry}}.star_idx'),
        'barcode_whitelist': get_barcode_whitelist(tech_type),
    }
    if tech_type == "10x_atac":
        input_['read'] = get_star_fastq_input(wc.cond, 'atac', 'read1')
        input_['mate'] = get_star_fastq_input(wc.cond, 'atac', 'read2')
        if config['preprocessing']['atac_rev_comp_barcode']:
            input_['barcode'] = get_star_fastq_input(wc.cond, 'atac', 'barcode_rc')
        else:
            input_['barcode'] = get_star_fastq_input(wc.cond, 'atac', 'barcode')
    else:
        input_['read'] = get_star_fastq_input(wc.cond, 'rna', 'read2')
        input_['barcode'] = get_star_fastq_input(wc.cond, 'rna', 'read1')
    return input_


def get_adapter_parameters(wc, input):
    tech_type = config['datasets'][wc.cond]['technology']
    whitelist = ' '.join(f'${{RELPATH}}/{fn}' for fn in input.barcode_whitelist)
    barcode_read_length = get_read_length(input.barcode[0])
    if tech_type == '10x_atac':
        params = f'''\
          --soloType "CB_samTagOut" \
          --soloCBwhitelist {whitelist} \
          --soloBarcodeReadLength 0 \
          --soloCBmatchWLtype "1MM" \
        '''
    elif tech_type.startswith('10x_rna'):
        params = f'''\
          --soloType "CB_UMI_Simple" \
          --soloCBwhitelist {whitelist} \
          --soloBarcodeReadLength {barcode_read_length} \
          --soloUMIdedup {config["alignment"]["star"]["rna_align_intron_max"]} \
          --soloCBmatchWLtype "1MM" \
          --soloCBlen 16 \
          --soloCBstart 1 \
          --soloUMIlen 12 \
          --soloUMIstart 17 \
        '''
    elif tech_type == 'bd_rna':
        params = f'''\
          --soloType "CB_UMI_Complex" \
          --soloCBwhitelist {whitelist} \
          --soloBarcodeReadLength {barcode_read_length} \
          --soloUMIdedup {config["alignment"]["star"]["rna_align_intron_max"]} \
          --soloAdapterSequence "NNNNNNNNNGTGANNNNNNNNNGACA" \
          --soloCBmatchWLtype "1MM" \
          --soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 \
          --soloUMIposition 3_10_3_17 \
        '''
    else:
        raise NotImplementedError()
    return params


def get_transform_flag(wc):
    dataset = config['datasets'][wc.cond]
    ref = dataset['reference_genotype']
    if wc.qry != ref:
        return '--genomeTransformOutput SAM'
    else:
        return ''


def get_temp_dir(wc, output):
    output_dir = os.path.split(output.bam)[0]
    tmp_dir = os.path.join(output_dir, f'{wc.cond}.{wc.qry}.tmpdir')
    return tmp_dir


def get_input_flags(wc, input):
    tech_type = config['datasets'][wc.cond]['technology']
    flag = '--readFilesCommand "zcat" --readFilesIn '
    flag += ','.join(f'${{RELPATH}}/{fn}' for fn in input.read)
    if tech_type == '10x_atac':
        flag += ' '
        flag += ','.join(f'${{RELPATH}}/{fn}' for fn in input.mate)
    flag += ' '
    flag += ','.join(f'${{RELPATH}}/{fn}' for fn in input.barcode)
    return flag
    

def get_spliced_alignment_params(wc):
    tech_type = config['datasets'][wc.cond]['technology']
    if tech_type == "10x_atac":
        return '''\
          --alignIntronMax 1 \
          --alignMatesGapMax {config["alignment"]["star"]["atac"]["mates_gap_max"]} \
        '''
    else:
        return f'''\
          --outFilterIntronMotifs RemoveNoncanonical \
          --alignSJoverhangMin 12 \
          --alignSJDBoverhangMin 4 \
          --alignIntronMin {config["alignment"]["star"]["rna_align_intron_min"]} \
          --alignIntronMax {config["alignment"]["star"]["rna_align_intron_max"]} \
        '''


rule STAR_consensus:
    input:
        unpack(STAR_consensus_input)
    output:
        bam=temp(results('aligned_data/haploid/{cond}.{qry}.sorted.bam')),
    params:
        star_tmp_dir=get_temp_dir,
        input_flag=get_input_flags,
        transform_flag=get_transform_flag,
        adapter_parameters=get_adapter_parameters,
        splicing_parameters=get_spliced_alignment_params,
        filter_multimap_nmax=config['alignment']['star']['filter_multimap_nmax'],
        filter_mismatch_nmax=config['alignment']['star']['filter_mismatch_nmax'],
        sort_mem=lambda wc, resources: (resources.mem_mb - 4096) * 1_000_000,
        n_files=lambda wc, threads: threads * 150 + 200,
    log:
        progress=results('logs/{cond}.{qry}.STAR_progress.log'),
        final=results('logs/{cond}.{qry}.STAR_final.log'),
        main=results('logs/{cond}.{qry}.STAR.log')
    threads: 24
    resources:
        mem_mb=lambda wildcards, threads: (threads + 4) * 2048,
    conda:
        '../env_yamls/star.yaml'
    shell:
        '''
        mkdir -p {params.star_tmp_dir}
        RELPATH=$(realpath --relative-to="{params.star_tmp_dir}" ".")
        cd {params.star_tmp_dir}
        ulimit -n {params.n_files}
        STAR \
          --runThreadN {threads} \
          --genomeDir "${{RELPATH}}/{input.index}" \
          {params.input_flag} \
          {params.adapter_parameters} \
          {params.splicing_parameters} \
          --outFilterMultimapNmax {params.filter_multimap_nmax} \
          --outFilterMismatchNmax {params.filter_mismatch_nmax} \
          --outSAMtype "BAM" "SortedByCoordinate" \
          --outBAMsortingBinsN 150 \
          --limitBAMsortRAM {params.sort_mem} \
          --outSAMattrRGline "ID:{wildcards.qry}" \
          --outSAMattributes "NH" "HI" "AS" "nM" "RG" "CB" "UB" \
          {params.transform_flag}

        cd $RELPATH
        mv {params.star_tmp_dir}/Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.star_tmp_dir}/Log.progress.out {log.progress}
        mv {params.star_tmp_dir}/Log.final.out {log.final}
        mv {params.star_tmp_dir}/Log.out {log.main}
          
        rm -rf {params.star_tmp_dir}
        '''


rule sort_bam_by_name:
    input:
        bam=results('aligned_data/haploid/{cond}.{qry}.sorted.bam')
    output:
        bam=temp(results('aligned_data/haploid/{cond}.{qry}.namesorted.bam'))
    threads: 12
    resources:
        mem_mb=20_000,
    conda:
        '../env_yamls/snco.yaml'
    shell:
        '''
        samtools sort -n -@ {threads} {input.bam} > {output.bam}
        '''


def get_merge_input(wc):
    dataset = config['datasets'][wc.cond]
    qrys = set()
    for geno in dataset['genotypes'].values():
        for qry in geno.values():
            qrys.add(qry)
    return {
        'bams': expand(results('aligned_data/haploid/{cond}.{qry}.namesorted.bam'),
                       cond=wc.cond, qry=qrys),
    }


rule merge_name_sorted_bams:
    input:
        unpack(get_merge_input)
    output:
        bam=temp(results('aligned_data/{cond}.namesorted.bam')),
    threads: 24
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
        bam=results('aligned_data/{cond}.namesorted.bam')
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
        collapse_ha_specific_alns.py \
          -o aligned_data/{wildcards.cond}.unsorted.bam {input.bam}
        samtools sort -@ {threads} \
          -T ${{TMPDIR}}/{wildcards.cond} \
          -o {output.bam} \
          aligned_data/{wildcards.cond}.unsorted.bam
        samtools index {output}
        rm aligned_data/{wildcards.cond}.unsorted.bam
        '''

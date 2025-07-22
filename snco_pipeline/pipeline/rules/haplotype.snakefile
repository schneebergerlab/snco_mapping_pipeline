def get_cb_tag_id(wc):
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type in ('takara_dna', 'plate_wgs'):
        return 'RG'
    else:
        return 'CB'


rule remove_pcr_duplicates:
    '''filter out PCR duplicates in single-cell ATAC libraries'''
    input:
        bam=results('aligned_data/{dataset_name}.sorted.bam'),
        bai=results('aligned_data/{dataset_name}.sorted.bam.bai'),
    output:
        bam=temp(results('aligned_data/{dataset_name}.deduped.bam')),
        bai=temp(results('aligned_data/{dataset_name}.deduped.bam.bai')),
    params:
        cb_tag=get_cb_tag_id,
    resources:
        mem_mb=10_000,
    conda:
        get_conda_env('htslib')
    threads: 12
    shell:
        format_command('''
        samtools markdup -r -@ {threads}
          --barcode-tag {params.cb_tag}
          {input.bam} {output.bam};

        samtools index {output.bam}
        ''')


def filter_input(wc):
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type not in ('10x_rna_v4', '10x_rna_v3', 'bd_rna') and \
            config['haplotyping']['preprocessing']['filter_pcr_duplicates']:
        input_file_type = 'deduped'
    else:
        # skips PCR deduplication by taking collapsed output as direct input
        input_file_type = 'sorted'
    return dict(
        bam=results(f'aligned_data/{wc.dataset_name}.{input_file_type}.bam'),
        bai=results(f'aligned_data/{wc.dataset_name}.{input_file_type}.bam.bai'),
    )


def get_filter_tag(wc):
    dataset = config['datasets'][wc.dataset_name]
    all_haplos = set()
    for geno in dataset['genotypes'].values():
        all_haplos.update(list(geno.values()))
    return ','.join(sorted(all_haplos))


rule filter_informative_reads:
    input:
        unpack(filter_input)
    output:
        bam=results('aligned_data/{dataset_name}.filtered.bam'),
        bai=results('aligned_data/{dataset_name}.filtered.bam.bai'),
        cb_whitelist=results('aligned_data/{dataset_name}.initial_whitelist.txt')
    conda:
        '../env_yamls/htslib.yaml'
    params:
        filter_tag=get_filter_tag,
        cb_tag=get_cb_tag_id,
        min_reads=config['haplotyping']['preprocessing']['min_informative_reads_per_barcode'],
    resources:
        mem_mb=20_000
    shell:
        format_command('''
        samtools view -b
          -e '[ha] != "{params.filter_tag}" && [{params.cb_tag}] != "-"'
          {input.bam}
        > {output.bam};

        samtools index {output.bam};

        samtools view --keep-tag "{params.cb_tag}" {output.bam} |
        awk '{{cb_count[substr($12, 6)]++}}
             END {{for (cb in cb_count)
                 {{if (cb_count[cb] > {params.min_reads})
                 {{print cb}}}}}}'
        > {output.cb_whitelist}
        ''')


def get_genotyping_params(wc):
    dataset = config['datasets'][wc.dataset_name]
    crosses = []
    for geno in dataset['genotypes'].values():
        crosses.append(':'.join(sorted(geno.values())))
    if len(crosses) > 1:
        genotype_params = '--genotype --clean-by-genotype -X '
    else:
        genotype_params = '-X '
    return genotype_params + ','.join(crosses)


def get_tech_specific_params(wc):
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type == "10x_atac":
        params = '''
          -x {x_flag}
          -y {ploidy}
          --cb-tag CB
          --hap-tag ha
          --hap-tag-type "multi_haplotype"
          --cb-correction-method "exact"
          --umi-collapse-method "none"
          --no-clean-bg
       '''
    elif tech_type in ("takara_dna", "plate_wgs"):
        params = '''
          -x {x_flag}
          -y {ploidy}
          --cb-tag RG
          --hap-tag ha
          --hap-tag-type "multi_haplotype"
          --cb-correction-method "none"
          --no-validate
          --umi-collapse-method "none"
          --no-clean-bg
          '''
        if tech_type == 'plate_wgs':
            params += '--no-predict-doublets'
    else:
        params = '''
          -x {x_flag}
          -y {ploidy}
          --cb-tag CB
          --umi-tag UB
          --hap-tag ha
          --hap-tag-type "multi_haplotype"
          --cb-correction-method "exact"
          --umi-collapse-method "exact"
       '''

    if tech_type.startswith('10x_rna'):
        x_flag = '10x_rna'
    elif tech_type == 'plate_wgs':
        x_flag = 'wgs'
    else:
        x_flag = tech_type
    params = params.format(
        x_flag=x_flag,
        ploidy=config['datasets'][wc.dataset_name]['ploidy']
    )

    return format_command(params.lstrip())


rule run_haplotyping:
    input:
        bam=results('aligned_data/{dataset_name}.filtered.bam'),
        bai=results('aligned_data/{dataset_name}.filtered.bam.bai'),
        cb_whitelist=results('aligned_data/{dataset_name}.initial_whitelist.txt')
    output:
        markers=results('haplotypes/{dataset_name}.markers.json'),
        preds=results('haplotypes/{dataset_name}.pred.json'),
        stats=results('haplotypes/{dataset_name}.pred.stats.tsv'),
    conda:
        get_conda_env('snco')
    params:
        bin_size=config['haplotyping']['snco']['genomic_bin_size'],
        rfactor=config['haplotyping']['snco']['segment_size'],
        term_rfactor=config['haplotyping']['snco']['terminal_segment_size'],
        cm_per_mb=config['haplotyping']['snco']['cm_per_mb'],
        genotyping_params=get_genotyping_params,
        tech_specific_params=get_tech_specific_params,
        min_reads_per_cb=config['haplotyping']['preprocessing']['min_informative_reads_per_barcode'],
        min_reads_per_chrom=config['haplotyping']['preprocessing']['min_informative_reads_per_chrom'],
        output_prefix=lambda wc: results(f'haplotypes/{wc.dataset_name}')
    threads: 32
    resources:
        mem_mb=100_000
    shell:
        format_command('''
        export OPENBLAS_NUM_THREADS=1;
        export OMP_NUM_THREADS=1;
        export MKL_NUM_THREADS=1;

        snco bam2pred -v debug -p {threads}
          --cb-whitelist-fn {input.cb_whitelist}
          -N {params.bin_size}
          -R {params.rfactor}
          -t {params.term_rfactor}
          -C {params.cm_per_mb}
          {params.tech_specific_params}
          --min-markers-per-cb {params.min_reads_per_cb}
          --min-markers-per-chrom {params.min_reads_per_chrom}
          --batch-size 128
          {params.genotyping_params}
          -o {params.output_prefix}
          {input.bam}
        ''')


rule haplotyping_report:
    input:
        markers=results('haplotypes/{dataset_name}.markers.json'),
        preds=results('haplotypes/{dataset_name}.pred.json'),
        stats=results('haplotypes/{dataset_name}.pred.stats.tsv'),
    log:
        notebook=results('analysis/{dataset_name}.haplotyping_report.py.ipynb')
    conda:
        get_conda_env('snco')
    notebook:
        '../notebook_templates/haplotyping_report.py.ipynb'
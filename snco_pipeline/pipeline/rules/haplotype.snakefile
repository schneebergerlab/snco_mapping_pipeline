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
        metrics=results('aligned_data/{dataset_name}.dedup_metrics.txt'),
        bai=temp(results('aligned_data/{dataset_name}.deduped.bam.bai')),
    params:
        cb_tag=get_cb_tag_id,
    resources:
        mem_mb=10_000,
    conda:
        '../env_yamls/picard.yaml'
    shell:
        '''
        samtools view -b -f2 {input.bam} > "{input.bam}.mapped.tmp.bam"
        picard MarkDuplicates \
          -I "{input.bam}.mapped.tmp.bam" \
          -M {output.metrics} \
          -O {output.bam} \
          --ASSUME_SORT_ORDER coordinate \
          --BARCODE_TAG "{params.cb_tag}" \
          --MAX_OPTICAL_DUPLICATE_SET_SIZE -1 \
          --REMOVE_DUPLICATES "true"
        rm "{input.bam}.mapped.tmp.bam"
        samtools index {output.bam}
        '''


def filter_input(wc):
    tech_type = config['datasets'][wc.dataset_name]['technology']
    if tech_type not in ('10x_rna_v4', '10x_rna_v3', 'bd_rna') and config['haplotyping']['filter_pcr_duplicates']:
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
        '''
        samtools view -b \
          -e '[ha] != "{params.filter_tag}" && [{params.cb_tag}] != "-"' \
          {input.bam} \
        > {output.bam}
        samtools index {output.bam}
        samtools view --keep-tag "{params.cb_tag}" {output.bam} | \
        awk '{{cb_count[substr($12, 6)]++}} \
             END {{for (cb in cb_count) \
                 {{if (cb_count[cb] > {params.min_reads}) \
                 {{print cb}}}}}}' \
        > {output.cb_whitelist}
        '''


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
    ploidy = config['datasets'][wc.dataset_name]['ploidy']
    if tech_type == "10x_atac":
        params = '''\
          -x 10x_atac \
          -y {ploidy} \
          --cb-tag CB \
          --hap-tag ha \
          --hap-tag-type "multi_haplotype" \
          --cb-correction-method "exact" \
          --umi-collapse-method "none" \
          --no-clean-bg \
       '''
    elif tech_type in ("takara_dna", "plate_wgs"):
        x_flag = "wgs" if tech_type == "plate_wgs" else tech_type
        params = f'''\
          -x {x_flag} \
          --no-validate \
          -y {ploidy} \
          --cb-tag RG \
          --hap-tag ha \
          --hap-tag-type "multi_haplotype" \
          --cb-correction-method "exact" \
          --umi-collapse-method "none" \
          --no-clean-bg \
        '''
        if tech_type == 'plate_wgs':
            params += '--no-predict-doublets'
    else:
        tech_type = '10x_rna' if tech_type.startswith('10x_rna') else 'bd_rna'
        params = f'''\
          -x {tech_type} \
          -y {ploidy} \
          --cb-tag CB \
          --umi-tag UB \
          --hap-tag ha \
          --hap-tag-type "multi_haplotype" \
          --cb-correction-method "exact" \
          --umi-collapse-method "exact" \
       '''
    return params


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
        '../env_yamls/snco.yaml'
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
        '''
        snco bam2pred -v debug -p {threads} \
          --cb-whitelist-fn {input.cb_whitelist} \
          -N {params.bin_size} \
          -R {params.rfactor} \
          -t {params.term_rfactor} \
          -C {params.cm_per_mb} \
          {params.tech_specific_params} \
          --min-markers-per-cb {params.min_reads_per_cb} \
          --min-markers-per-chrom {params.min_reads_per_chrom} \
          --batch-size 128 \
          {params.genotyping_params} \
          -o {params.output_prefix} \
          {input.bam}
        '''

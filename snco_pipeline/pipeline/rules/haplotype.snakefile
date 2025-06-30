rule remove_pcr_duplicates:
    input:
        bam=results('aligned_data/{cond}.sorted.bam'),
        bai=results('aligned_data/{cond}.sorted.bam.bai'),
    output:
        bam=temp(results('aligned_data/{cond}.deduped.bam')),
        metrics=results('aligned_data/{cond}.dedup_metrics.txt'),
        bai=temp(results('aligned_data/{cond}.deduped.bam.bai')),
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
          --BARCODE_TAG "CB" \
          --MAX_OPTICAL_DUPLICATE_SET_SIZE -1 \
          --REMOVE_DUPLICATES "true"
        rm "{input.bam}.mapped.tmp.bam"
        samtools index {output.bam}
        '''


def filter_input(wc):
    tech_type = config['datasets'][wc.cond]['technology']
    input_file_type = 'deduped' if tech_type == '10x_atac' else 'sorted'
    return dict(
        bam=results(f'aligned_data/{wc.cond}.{input_file_type}.bam'),
        bai=results(f'aligned_data/{wc.cond}.{input_file_type}.bam.bai'),
    )


def get_filter_tag(wc):
    dataset = config['datasets'][wc.cond]
    all_haplos = set()
    for geno in dataset['genotypes'].values():
        all_haplos.update(list(geno.values()))
    return ','.join(sorted(all_haplos))


rule filter_informative_reads:
    input:
        unpack(filter_input)
    output:
        bam=results('aligned_data/{cond}.filtered.bam'),
        bai=results('aligned_data/{cond}.filtered.bam.bai'),
        cb_whitelist=results('aligned_data/{cond}.initial_whitelist.txt')
    conda:
        '../env_yamls/pysam.yaml'
    params:
        filter_tag=get_filter_tag,
        min_informative_reads_per_cb=100,
        bc_len=lambda wc: 16 if config['datasets'][wc.cond]['technology'] != 'bd' else 29
    resources:
        mem_mb=20_000
    shell:
        '''
        samtools view -b \
          -e '[ha] != "{params.filter_tag}" && [CB] != "-"' \
          {input.bam} \
        > {output.bam}
        samtools index {output.bam}
        samtools view --keep-tag "CB" {output.bam} | \
        awk '{{cb_count[substr($12, 6, {params.bc_len})]++}} \
             END {{for (cb in cb_count) \
                 {{if (cb_count[cb] > {params.min_informative_reads_per_cb}) \
                 {{print cb}}}}}}' \
        > {output.cb_whitelist}
        '''


def get_genotyping_params(wc):
    dataset = config['datasets'][wc.cond]
    crosses = []
    for geno in dataset['genotypes'].values():
        crosses.append(':'.join(sorted(geno.values())))
    if len(crosses) > 1:
        genotype_params = '--genotype --clean-by-genotype -X '
    else:
        genotype_params = '-X '
    return genotype_params + ','.join(crosses)


def get_tech_specific_params(wc):
    tech_type = config['datasets'][wc.cond]['technology']
    if tech_type == "10x_atac":
        params = '''
          -x 10x_atac \
          -y haploid \
          --cb-tag CB \
          --hap-tag ha \
          --hap-tag-type "multi_haplotype" \
          --cb-correction-method "exact" \
          --umi-collapse-method "none" \
          --no-clean-bg \
       '''
    else:
        tech_type == '10x_rna' if tech_type.startswith('10x_rna') else 'bd_rna'
        params = f'''
          -x {tech_type} \
          -y haploid \
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
        bam=results('aligned_data/{cond}.filtered.bam'),
        bai=results('aligned_data/{cond}.filtered.bam.bai'),
        cb_whitelist=results('aligned_data/{cond}.initial_whitelist.txt')
    output:
        markers=results('haplotypes/{cond}.markers.json'),
        preds=results('haplotypes/{cond}.pred.json'),
        stats=results('haplotypes/{cond}.pred.stats.tsv'),
    conda:
        '../env_yamls/snco.yaml'
    params:
        bin_size=25_000,
        genotyping_params=get_genotyping_params,
        tech_specific_params=get_tech_specific_params,
    threads: 32
    resources:
        mem_mb=100_000
    shell:
        '''
        snco bam2pred -v debug -p {threads} \
          --cb-whitelist-fn {input.cb_whitelist} \
          -N {params.bin_size} \
          {params.tech_specific_params} \
          --min-markers-per-cb 100 \
          --min-markers-per-chrom 20 \
          --batch-size 128 \
          {params.genotyping_params} \
          -o "haplotypes/{wildcards.cond}" \
          {input.bam}
        '''

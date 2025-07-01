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


def get_cb_tag_id(wc):
    tech_type = config['datasets'][wc.cond]['technology']
    if tech_type in ('takara_dna', 'plate_wgs'):
        return 'RG'
    else:
        return 'CB'

rule filter_informative_reads:
    input:
        unpack(filter_input)
    output:
        bam=results('aligned_data/{cond}.filtered.bam'),
        bai=results('aligned_data/{cond}.filtered.bam.bai'),
        cb_whitelist=results('aligned_data/{cond}.initial_whitelist.txt')
    conda:
        '../env_yamls/bcftools.yaml'
    params:
        filter_tag=get_filter_tag,
        cb_tag=get_cb_tag_id,
        min_informative_reads_per_cb=config['haplotyping']['min_informative_reads_per_barcode'],
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
    ploidy = config['datasets'][wc.cond]['ploidy']
    if tech_type == "10x_atac":
        params = '''
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
        tech_type = "wgs" if tech_type == "plate_wgs" else tech_type
        params = f'''
          -x {tech_type} \
          -y {ploidy} \
          --cb-tag RG \
          --hap-tag ha \
          --hap-tag-type "multi_haplotype" \
          --cb-correction-method "exact" \
          --umi-collapse-method "none" \
          --no-clean-bg \
       '''
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
        min_reads_per_cb=config['haplotyping']['min_informative_reads_per_barcode'],
        min_reads_per_chrom=config['haplotyping']['min_informative_reads_per_chrom'],
    threads: 32
    resources:
        mem_mb=100_000
    shell:
        '''
        snco bam2pred -v debug -p {threads} \
          --cb-whitelist-fn {input.cb_whitelist} \
          -N {params.bin_size} \
          {params.tech_specific_params} \
          --min-markers-per-cb {params.min_reads_per_cb} \
          --min-markers-per-chrom {params.min_reads_per_chrom} \
          --batch-size 128 \
          {params.genotyping_params} \
          -o "haplotypes/{wildcards.cond}" \
          {input.bam}
        '''

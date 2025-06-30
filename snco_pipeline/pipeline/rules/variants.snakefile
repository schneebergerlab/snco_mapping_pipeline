def get_wga_input(wc):
    return {
        'ref': get_fasta(wc.ref),
        'qry': get_fasta(wc.qry),
    }


rule minimap2_wga:
    input:
        unpack(get_wga_input)
    output:
        annotation('wga/{ref}.{qry}.bam')
    threads:
        12
    resources:
        mem_mb=lambda wc, threads: 1024 * threads,
    conda:
        '../env_yamls/minimap2.yaml'
    shell:
        '''
        minimap2 --eqx -t 8 -ax asm10 -z500 {input.ref} {input.qry} | \
        samtools view -bS | \
        samtools sort -o - - > {output}
        samtools index {output}
        '''


def get_syri_input(wc):
    input_ = get_wga_input(wc)
    input_['bam'] = annotation('wga/{ref}.{qry}.bam')
    return input_


rule run_syri:
    input:
        unpack(get_syri_input)
    output:
        vcf=annotation('vcf/syri/{ref}.{qry}.syri.vcf'),
        syri=annotation('vcf/syri/{ref}.{qry}.syri.out'),
    resources:
        mem_mb=lambda wc, threads: 2048 * threads,
    conda:
        '../env_yamls/msyd.yaml'
    shell:
        '''
        syri -F B -f --hdrseq --dir annotation/vcf/syri --prefix {wildcards.ref}.{wildcards.qry}. \
          -c {input.bam} -q {input.qry} -r {input.ref} --samplename {wildcards.qry}
        '''


rule filter_syri_snps_for_star_consensus:
    input:
        vcf=annotation('vcf/syri/{ref}.{qry}.syri.vcf')
    output:
        vcf=annotation('vcf/star_consensus/syri/{ref}.{qry}.vcf')
    conda:
        '../env_yamls/snco.yaml'
    resources:
        mem_mb=30_000
    shell:
        '''
        syri_vcf_to_stardiploid.py -M0 -i50 -n {wildcards.qry} {input} {output}
        '''


def get_msyd_input(wc):
    ref, *qry_names = wc.geno_group.split('_')
    return {
        'bams': expand(annotation('wga/{ref}.{qry}.bam'),
                       ref=ref, qry=qry_names),
        'syri': expand(annotation('vcf/syri/{ref}.{qry}.syri.out'),
                       ref=ref, qry=qry_names),
        'vcf': expand(annotation('vcf/syri/{ref}.{qry}.syri.vcf'),
                       ref=ref, qry=qry_names),
        'qry_fasta': [get_fasta(qry) for qry in qry_names],
        'ref_fasta': get_fasta(ref)
    }


rule msyd_input:
    input:    
        unpack(get_msyd_input)
    output:
        cfg=temp(annotation('vcf/msyd/{geno_group}.msyd_config.tsv'))
    run:
        ref, *qry_names = wildcards.geno_group.split('_')
        with open(output.cfg, 'w') as f:
            f.write('#name\taln\tsyri\tvcf\tgenome\n')
            for qry in qry_names:
                f.write(
                    f'{qry}\t'
                    f'annotation/wga/{ref}.{qry}.bam\t'
                    f'annotation/vcf/syri/{ref}.{qry}.syri.out\t'
                    f'annotation/vcf/syri/{ref}.{qry}.syri.vcf\t'
                    f'{get_genome_fasta(qry)}\n'
                )


rule run_msyd:
    input:
        unpack(get_msyd_input),
        cfg=annotation('vcf/msyd/{geno_group}.msyd_config.tsv'),
    output:
        pff=annotation('vcf/msyd/{geno_group,\w+}.pansyn.pff'),
        vcf=annotation('vcf/msyd/{geno_group,\w+}.vcf'),
    conda:
        '../env_yamls/msyd.yaml'
    shell:
        '''
        msyd call --core -i {input.cfg} -r {input.ref_fasta} -o {output.pff} -m {output.vcf}
        '''


rule filter_msyd_snps_for_star_consensus:
    input:
        vcf=annotation('vcf/msyd/{geno_group}.vcf')
    output:
        vcf=annotation('vcf/star_consensus/msyd/{geno_group}.{qry}.vcf'),
    conda:
        '../env_yamls/pysam.yaml'
    params:
        max_indel_size=config['variants']['max_indel_size'],
    shell:
        r'''
        bcftools annotate \
          --exclude 'ALT ~ "CORESYN"' \
          --remove "FORMAT/CHR,FORMAT/START,FORMAT/END,INFO/PID" \
          {input.vcf} | \
        grep -v "^##ALT" | grep -v "^##INFO" | \
        bcftools sort | \
        bcftools filter -S0 -e 'GT=="."' | \
        bcftools view -s {wildcards.qry} | \
        bcftools view -i 'GT=="alt"' | \
        bcftools view -G \
          -e "STRLEN(REF)>{params.max_indel_size} || STRLEN(ALT)>{params.max_indel_size}" \
        > {output.vcf}
        '''

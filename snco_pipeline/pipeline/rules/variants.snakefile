def get_wga_input(wc):
    '''input genomes for minimap2 and syri'''
    return {
        'ref': get_fasta(wc.ref),
        'qry': get_fasta(wc.qry),
    }


rule minimap2_wga:
    '''align two chromosome-scale genome assemblies with minimap2 for analysis with syri'''
    input:
        unpack(get_wga_input)
    output:
        annotation('wga/{ref}.{qry}.bam')
    threads:
        12
    resources:
        mem_mb=lambda wc, threads: 1024 * threads,
    params:
        preset=config['variants']['minimap2']['preset'],
        zdrop=config['variants']['minimap2']['zdrop'],
    conda:
        '../env_yamls/minimap2.yaml'
    shell:
        '''
        minimap2 --eqx -t {threads} \
          -ax {params.preset} \
          -z{params.zdrop} \
          {input.ref} {input.qry} | \
        samtools view -bS | \
        samtools sort -o - - > {output}
        samtools index {output}
        '''


def get_syri_input(wc):
    '''full syri input, requires the two genomes plus their minimap2 alignment (bam file)'''
    input_ = get_wga_input(wc)
    input_['bam'] = annotation('wga/{ref}.{qry}.bam')
    return input_


rule run_syri:
    '''use syri to convert a wga of two genome assemblies into a  vcf file of synteny, as well as snps/indels'''
    input:
        unpack(get_syri_input)
    output:
        vcf=annotation('vcf/syri/{ref}.{qry}.syri.vcf'),
        syri=annotation('vcf/syri/{ref}.{qry}.syri.out'),
    resources:
        mem_mb=lambda wc, threads: 1024 * threads,
    conda:
        '../env_yamls/msyd.yaml'
    params:
        filter_alns='' if config['variants']['syri']['use_low_qual_filters'] else '-f',
        out_dir=annotation('vcf/syri')
    shell:
        '''
        syri -F B --hdrseq \
          {params.filter_alns} \
          --dir {params.out_dir} \
          --prefix {wildcards.ref}.{wildcards.qry}. \
          -c {input.bam} \
          -q {input.qry} \
          -r {input.ref} \
          --samplename {wildcards.qry}
        '''


def get_msyd_input(wc):
    '''full input for msyd, requires genomes, wga bams, and syri output as vcf/syri.out files'''
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
    '''builds a config file for msyd'''
    input:    
        unpack(get_msyd_input)
    output:
        cfg=temp(annotation('vcf/msyd/{geno_group}.msyd_config.tsv'))
    params:
        annot=config['annotation_dir']
    run:
        ref, *qry_names = wildcards.geno_group.split('_')
        with open(output.cfg, 'w') as f:
            f.write('#name\taln\tsyri\tvcf\tgenome\n')
            for qry in qry_names:
                f.write(
                    f'{qry}\t'
                    f'{params.annot}/wga/{ref}.{qry}.bam\t'
                    f'{params.annot}/vcf/syri/{ref}.{qry}.syri.out\t'
                    f'{params.annot}/vcf/syri/{ref}.{qry}.syri.vcf\t'
                    f'{get_fasta(qry)}\n'
                )


rule run_msyd:
    '''
    runs msyd to identify "core" synteny between two or more assemblies, 
    outputs a vcf of SNPs/indels in core regions
    '''
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
        msyd call --core \
          -i {input.cfg} \
          -r {input.ref_fasta} \
          -o {output.pff} \
          -m {output.vcf}.tmp.vcf
        bcftools annotate \
          --exclude 'ALT ~ "CORESYN"' \
          --remove "FORMAT/CHR,FORMAT/START,FORMAT/END,INFO/PID" \
          {output.vcf}.tmp.vcf | \
        grep -v "^##ALT" | grep -v "^##INFO" | \
        bcftools sort | \
        bcftools filter -S0 -e 'GT=="."' > {output.vcf}
        rm {output.vcf}.tmp.vcf
        '''


def filter_snps_vcf_input(wc):
    '''input for vcf filtering, can be either msyd output created from wga or a user defined vcf file'''
    vcf_fns = config['annotations']['vcf_fns']
    ref, *qry_names = wc.geno_group.split('_')
    if ref in vcf_fns and vcf_fns[ref] is not None:
        return ancient(annotation(vcf_fns[ref]))
    else:
        return annotation('vcf/msyd/{geno_group}.vcf')


rule filter_msyd_snps_for_star_consensus:
    '''
    filter and reformat a vcf file (msyd output or user defined) to produce a VCF
    for a single query vs reference comparison suitable for STAR consensus
    '''
    input:
        vcf=filter_snps_vcf_input
    output:
        vcf=annotation('vcf/star_consensus/msyd/{geno_group}.{qry}.vcf'),
    conda:
        '../env_yamls/htslib.yaml'
    params:
        max_indel_size=config['variants']['star_consensus']['max_indel_size'],
    shell:
        r'''
        bcftools view -s {wildcards.qry} {input.vcf} | \
        bcftools view -i 'GT=="alt"' | \
        bcftools view -G \
          -e "STRLEN(REF)>{params.max_indel_size} || STRLEN(ALT)>{params.max_indel_size}" \
        > {output.vcf}
        '''


ORG_MAP = {'homo_sapiens': 'hsa', 'mus_musculus': 'mmu', 'rattus_norvegicus': 'rnor'}

MIRTRACE_DIR = join(

rule mirtrace_config:
    input:
        fastq = expand(get_filtered_fastq, sample=SAMPLES)    
    output:
        

rule mirtrace_qc:
    input:
        rules.mirtrace_config.output
    params:
        org = ORG_MAP[config['organism']]
    singularity:
    
    shell:
        'mirtrace qc '
        '--config {input.config} '
        '--species {params.org} '
        '--adapter {params.adapter} '
        '--protocol {params.protocol} '
        '--output-dir {params.outdir} '
        '--force '
        
    

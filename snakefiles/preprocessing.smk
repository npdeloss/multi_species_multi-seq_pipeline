# Import naming conventions for reads
configfile: 'configs/reads.yaml'
# Import preprocessing default configuration
configfile: 'configs/preprocessing.yaml'

# Import user configuration
configfile: 'config.yaml'

# Establish ruleorder: paired if available, else single end
# ruleorder: preprocessing_fastp_paired_end > preprocessing_fastp_single_end
# ruleorder: preprocessing_fastp_paired_end > preprocessing_bbmap_single_end
ruleorder: preprocessing_fastp_paired_end > preprocessing_homertools_single_end

# Import naming conventions for reads as variables
reads_prefix = config['reads']['prefix']
reads1_suffix = config['reads']['reads1_suffix']
reads2_suffix = config['reads']['reads2_suffix']
fastq_suffix = config['reads']['fastq_suffix']
compression_suffix = config['reads']['compression_suffix']

# Construct options string for paired end processing
preprocessing_fastp_prefix = config['preprocessing']['fastp']['prefix']
preprocessing_fastp_paired_end_options = config['preprocessing']['fastp']['paired_end']['options']
adapter_r1 = config['preprocessing']['fastp']['paired_end']['adapter_r1']
adapter_r1_option = f'--adapter_sequence {adapter_r1}' if adapter_r1 != '' else ''
adapter_r2 = config['preprocessing']['fastp']['paired_end']['adapter_r1']
adapter_r2_option = f'--adapter_sequence_r2 {adapter_r2}' if adapter_r2 != '' else ''
preprocessing_fastp_paired_end_options = f'{preprocessing_fastp_paired_end_options} {adapter_r1_option} {adapter_r2_option}'

# Construct options string for single end processing
preprocessing_fastp_single_end_options = config['preprocessing']['fastp']['single_end']['options']
adapter = config['preprocessing']['fastp']['single_end']['adapter']
adapter_option = f'-a {adapter}' if adapter != '' else ''
preprocessing_fastp_single_end_options = f'{preprocessing_fastp_single_end_options} {adapter_option}'

# Import number of threads

preprocessing_fastp_threads = config['preprocessing']['fastp']['threads']

# Paired end
rule preprocessing_fastp_paired_end:
    input:
        reads1 = reads_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        reads2 = reads_prefix + '{basename}' + reads2_suffix + fastq_suffix + compression_suffix
    output:
        reads1 = preprocessing_fastp_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        reads2 = preprocessing_fastp_prefix + '{basename}' + reads2_suffix + fastq_suffix + compression_suffix,
        html = preprocessing_fastp_prefix + '{basename}' + '.html',
        json = preprocessing_fastp_prefix + '{basename}' + '.json'
    log:
        preprocessing_fastp_prefix + '{basename}' + '.log'
    threads:
        preprocessing_fastp_threads
    params:
        options = preprocessing_fastp_paired_end_options
    conda:
        '../envs/fastp.yaml'
    shell:
        """
        mkdir -p $(dirname {output.html})
        fastp \
        -i {input.reads1} \
        -I {input.reads2} \
        -o {output.reads1} \
        -O {output.reads2} \
        -h {output.html} \
        -j {output.json} \
        -z 9 \
        {params.options} \
        -w {threads} \
        &> {log}
        """

# Single end
# rule preprocessing_fastp_single_end:
#     input:
#         reads1 = reads_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix
#     output:
#         reads1 = preprocessing_fastp_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
#         html = preprocessing_fastp_prefix + '{basename}' + '.html',
#         json = preprocessing_fastp_prefix + '{basename}' + '.json'
#     log:
#         preprocessing_fastp_prefix + '{basename}' + '.log'
#     threads:
#         preprocessing_fastp_threads
#     params:
#         options = preprocessing_fastp_single_end_options
#     conda:
#         '../envs/fastp.yaml'
#     shell:
#         """
#         mkdir -p $(dirname {output.html})
#         fastp \
#         -i {input.reads1} \
#         -o {output.reads1} \
#         -h {output.html} \
#         -j {output.json} \
#         -z 9 \
#         {params.options} \
#         -w {threads} \
#         &> {log}
#         """

# rule preprocessing_bbmap_single_end:
#     input:
#         reads1 = reads_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix, 
#         ref = 'resources/adapters.fa'
#     output:
#         reads1 = preprocessing_fastp_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix
#     log:
#         preprocessing_fastp_prefix + '{basename}' + '.log'
#     threads:
#         preprocessing_fastp_threads
#     conda:
#         '../envs/bbmap.yaml'
#     shell:
#         """
#         mkdir -p $(dirname {output.reads1})
#         bbduk.sh -Xmx10g \
#         in={input.reads1} \
#         out={output.reads1} \
#         ref={input.ref} \
#         threads={threads} \
#         ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20 \
#         &> {log}
#         """

rule preprocessing_homertools_single_end:
    input:
        reads1 = reads_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix
    output:
        reads1 = preprocessing_fastp_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output.reads1})
        cp {input.reads1} {output.reads1}
        homerTools trim -3 AGATCGGAAGAGCACACGTCT -mis 2 -minMatchLength 4 -min 20 {output.reads1}
        cat {output.reads1}.trimmed | gzip -c > {output.reads1}
        """

# Import naming conventions for reads as variables
reads1_suffix = config['reads']['reads1_suffix']
reads2_suffix = config['reads']['reads2_suffix']
fastq_suffix = config['reads']['fastq_suffix']
compression_suffix = config['reads']['compression_suffix']

# Import prefixes for input and output
qc_fastqc_prefix = config['qc']['fastqc']['prefix']
qc_fastqc_input_prefix = config['qc']['fastqc']['input_prefix']

rule qc_fastqc:
    input:
        reads = qc_fastqc_input_prefix + '{basename_and_read_suffix}' + fastq_suffix + compression_suffix
    output:
        html = qc_fastqc_prefix + '{basename_and_read_suffix}_fastqc.html',
        zip = qc_fastqc_prefix + '{basename_and_read_suffix}_fastqc.zip'
    log:
        qc_fastqc_prefix + '{basename_and_read_suffix}_fastqc.log'
    threads:
        1
    conda:
        '../envs/fastqc.yaml'
    shell:
        """
        mkdir -p $(dirname {output.html})
        fastqc -o $(dirname {output.html}) {input.reads} &> {log}
        """

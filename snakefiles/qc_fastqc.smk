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

rule qc_fastqc_unzip:
    input:
        qc_fastqc_prefix + '{basename_and_read_suffix}_fastqc.zip'
    output:
        qc_fastqc_prefix + '{basename_and_read_suffix}_fastqc/summary.txt'
    log:
        qc_fastqc_prefix + '{basename_and_read_suffix}_fastqc.unzip.log'
    shell:
        """
        cd $(dirname {input})
        unzip $(basename {input}) &> $(basename {log})
        """

def get_all_basenames_and_read_suffixes():
    reads1_basenames = get_read_basenames_for_organism(None, 
                                                       library_type = None, 
                                                       basename_pattern = basename_pattern, 
                                                       reads_prefix = reads_prefix, 
                                                       suffix = reads1_suffix + fastq_suffix + compression_suffix)
    reads2_basenames = get_read_basenames_for_organism(None, 
                                                   library_type = None, 
                                                   basename_pattern = basename_pattern, 
                                                   reads_prefix = reads_prefix, 
                                                   suffix = reads1_suffix + fastq_suffix + compression_suffix)
    basenames_and_read_suffixes = [basename + reads1_suffix for basename in reads1_basenames]
    basenames_and_read_suffixes = basenames_and_read_suffixes + [basename + reads2_suffix for basename in reads2_basenames]
    return basenames_and_read_suffixes

rule qc_fastqc_index_from_reads:
    input:
        [qc_fastqc_prefix + basename_and_read_suffix + '/summary.txt' for basename_and_read_suffix in get_all_basenames_and_read_suffixes()]
    output:
        qc_fastqc_prefix + 'index.txt'
    run:
        summary_list = [summary for summary in input]
        qc_dirs_list = sorted(['/'.join(summary.split('/')[:-1])+'/' for summary in summary_list])
        qc_dirs_string = '\n'.join(qc_dirs_list)+'\n'
        with (output[0], 'w') as outfile:
            outfile.write(qc_dirs_string)
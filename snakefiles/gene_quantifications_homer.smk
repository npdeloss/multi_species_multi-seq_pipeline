gq_homer = config['gene_quantifications']['homer']
gq_homer_prefix = gq_homer['prefix']
gq_homer_input_prefix = gq_homer['input_prefix']
gq_homer_threads = gq_homer['threads']

gq_homer_index = config['gene_quantifications']['homer']['index']
gq_homer_index_input_prefix = gq_homer_index['input_prefix']

gq_homer_paired_end_options = gq_homer['paired_end']['params']['options']
gq_homer_single_end_options = gq_homer['single_end']['params']['options']

ruleorder: gene_quantifications_homer_paired_end > gene_quantifications_homer_single_end:

rule gene_quantifications_homer_paired_end:
    input:
        tag_dir = gq_kallisto_input_prefix + '{basename}' + '/tagInfo.txt',
        reads2 = reads_prefix + '{basename}' + reads2_suffix + fastq_suffix + compression_suffix,
        annotation_gtf = gq_homer_index_input_prefix + '/annotation.gtf'
    output:
        gq_homer_prefix + '{basename}' + '/{norm_method}.txt'
    log:
        gq_homer_prefix + '{basename}' + '/{norm_method}.log'
    params:
        options = gq_homer_paired_end_options
    shadow:
        'shallow'
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        analyzeRepeats.pl {input.annotation_gtf} none {params.options} -d $(dirname {input.tag_dir}) &> {log}
        """

rule gene_quantifications_homer_single_end:
    input:
        tag_dir = gq_kallisto_input_prefix + '{basename}' + '/tagInfo.txt',
        annotation_gtf = gq_homer_index_input_prefix + '/annotation.gtf'
    output:
        gq_homer_prefix + '{basename}' + '/{norm_method}.txt'
    log:
        gq_homer_prefix + '{basename}' + '/{norm_method}.log'
    params:
        options = gq_homer_single_end_options
    shadow:
        'shallow'
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        analyzeRepeats.pl {input.annotation_gtf} none {params.options} -d $(dirname {input.tag_dir}) &> {log}
        """

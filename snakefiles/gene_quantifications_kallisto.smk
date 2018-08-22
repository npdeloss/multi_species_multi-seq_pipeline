gq_kallisto = config['gene_quantifications']['kallisto']
gq_kallisto_prefix = gq_kallisto['prefix']
gq_kallisto_input_prefix = gq_kallisto['input_prefix']
gq_kallisto_threads = gq_kallisto['threads']

gq_kallisto_index = config['gene_quantifications']['kallisto']['index']
gq_kallisto_index_prefix = gq_kallisto_index['prefix']
gq_kallisto_index_options = gq_kallisto_index['params']['options']
gq_kallisto_index_input_prefix = gq_kallisto_index['input_prefix']

gq_kallisto_paired_end_options = gq_kallisto['paired_end']['params']['options']
gq_kallisto_single_end_options = gq_kallisto['paired_end']['params']['options']

rule index_kallisto:
    input:
        gq_kallisto_index_input_prefix + 'transcriptome.fa'
    output:
        gq_kallisto_index_prefix + 'transcriptome.idx'
    log:
        gq_kallisto_index_prefix + 'transcriptome.log'
    params:
        options = gq_kallisto_index_options
    conda:
        '../envs/kallisto.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        kallisto index {params.options} -i {output} {input} &> {log}
        """

ruleorder: quantifications_kallisto_paired_end > quantifications_kallisto_single_end

rule gene_quantifications_kallisto_paired_end:
    input:
        reads1 = gq_kallisto_input_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        reads2 = gq_kallisto_input_prefix + '{basename}' + reads2_suffix + fastq_suffix + compression_suffix,
        index = gq_kallisto_index_prefix + 'transcriptome.idx'
    output:
        abundance_tsv = gq_kallisto_prefix + '{basename}/abundance.tsv',
        abundance_h5 = gq_kallisto_prefix + '{basename}/abundance.h5',
        runinfo = gq_kallisto_prefix + '{basename}/run_info.json'
    log:
        gq_kallisto_prefix + '{basename}.log'
    threads:
        gq_kallisto_threads
    params:
        options = gq_kallisto_paired_end_options
    shell:
        """
        mkdir -p $(dirname {output.abundance_tsv})
        kallisto quant \
        -i {input.index} \
        -o $(dirname {output.abundance_tsv}) \
        {params.options} \
        {input.reads1} \
        {input.reads2}
        """

rule gene_quantifications_kallisto_single_end:
    input:
        reads1 = gq_kallisto_input_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        index = gq_kallisto_index_prefix + 'transcriptome.idx'
    output:
        abundance_tsv = gq_kallisto_prefix + '{basename}/abundance.tsv',
        abundance_h5 = gq_kallisto_prefix + '{basename}/abundance.h5',
        runinfo = gq_kallisto_prefix + '{basename}/run_info.json'
    log:
        gq_kallisto_prefix + '{basename}.log'
    threads:
        gq_kallisto_threads
    params:
        options = gq_kallisto_single_end_options
    shell:
        """
        mkdir -p $(dirname {output.abundance_tsv})
        kallisto quant \
        -i {input.index} \
        -o $(dirname {output.abundance_tsv}) \
        {params.options} \
        --single {input.reads1}
        """

rule gene_quantifications_kallisto_library_type_index:
    input:
        lambda wildcards: [(gq_kallisto_prefix + basename + '/abundance.tsv').format(**wildcards) for basename in get_read_basenames_for_organism(wildcards.organism, 
                                                                                                                                                        library_type = wildcards.library_type)]
    output:
        gq_kallisto_prefix + 'index.{library_type}.txt'
    run:
        shell(f'mkdir -p $(dirname {output[0]})')
        input_files = [input_file for input_file in input]
        input_dirs = sorted(['/'.join(input_file.split('/')[:-1])+'/' for input_file in input_files])
        input_dirs_string = '\n'.join(input_dirs)+'\n'
        with open(output[0], 'w') as outfile:
            outfile.write(input_dirs_string)

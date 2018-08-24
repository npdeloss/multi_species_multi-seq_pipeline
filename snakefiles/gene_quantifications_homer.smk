gq_homer = config['gene_quantifications']['homer']
gq_homer_prefix = gq_homer['prefix']
gq_homer_input_prefix = gq_homer['input_prefix']

gq_homer_index = config['gene_quantifications']['homer']['index']
gq_homer_index_input_prefix = gq_homer_index['input_prefix']

gq_homer_paired_end_options = gq_homer['paired_end']['params']['options']
gq_homer_single_end_options = gq_homer['single_end']['params']['options']

def get_gw_homer_option_for_norm_method(norm_method):
    if norm_method.startswith('norm'):
        split_norm_method = norm_method.split('_')
        if len(split_norm_method) is 2:
            return '-'+' '.join(split_norm_method)
        else:
            return '-norm'
    elif norm_method in ['raw', 'rpkm', 'log', 'sqrt']:
        return '-'+norm_method
    else:
        return '-raw'

ruleorder: gene_quantifications_homer_paired_end > gene_quantifications_homer_single_end

rule gene_quantifications_homer_paired_end:
    input:
        tag_dir = gq_homer_input_prefix + '{basename}' + '/tagInfo.txt',
        reads2 = reads_prefix + '{basename}' + reads2_suffix + fastq_suffix + compression_suffix,
        annotation_gtf = gq_homer_index_input_prefix + 'annotation.gtf'
    output:
        gq_homer_prefix + '{basename}' + '/{norm_method}.txt'
    log:
        gq_homer_prefix + '{basename}' + '/{norm_method}.log'
    params:
        options = gq_homer_paired_end_options
        norm_option = lambda wildcards: get_gw_homer_options_for_norm_method(wildcards.norm_method)
    shadow:
        'shallow'
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        analyzeRepeats.pl \
        {input.annotation_gtf} \
        none \
        {params.options} \
        {params.norm_option} \
        -d $(dirname {input.tag_dir}) \
        > {output} \
        2>> {log}
        """

rule gene_quantifications_homer_single_end:
    input:
        tag_dir = gq_kallisto_input_prefix + '{basename}' + '/tagInfo.txt',
        annotation_gtf = gq_homer_index_input_prefix + 'annotation.gtf'
    output:
        gq_homer_prefix + '{basename}' + '/{norm_method}.txt'
    log:
        gq_homer_prefix + '{basename}' + '/{norm_method}.log'
    params:
        options = gq_homer_single_end_options
        norm_option = lambda wildcards: get_gw_homer_options_for_norm_method(wildcards.norm_method)
    shadow:
        'shallow'
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        analyzeRepeats.pl \
        {input.annotation_gtf} \
        none \
        {params.options} \
        {params.norm_option} \
        -d $(dirname {input.tag_dir}) \
        > {output} \
        2>> {log}
        """

rule gene_quantifications_homer:
    input:
        gq_homer_prefix + '{basename}' + '/{norm_method}.txt'
    output:
        gq_homer_prefix + '{basename}' + '/{norm_method}.tsv'
    shell:
        """
        echo -e 'gene_id\t{wildcards.norm_method}' > {output}
        cat {input} | cut -f1,9 | tail -n +2 >> {output}
        """

rule gene_quantifications_homer_library_type_index:
    input:
        lambda wildcards: [(gq_homer_prefix + basename + '/{norm_method}.tsv').format(**wildcards) for basename in get_read_basenames_for_organism(wildcards.organism, 
                                                                                                                                                        library_type = wildcards.library_type)]
    output:
        gq_homer_prefix + 'index.{norm_method}.{library_type}.txt'
    run:
        shell(f'mkdir -p $(dirname {output[0]})')
        input_files = [input_file for input_file in input]
        input_dirs = sorted(['/'.join(input_file.split('/')[:-1])+'/' for input_file in input_files])
        input_dirs_string = '\n'.join(input_dirs)+'\n'
        with open(output[0], 'w') as outfile:
            outfile.write(input_dirs_string)

rule gene_quantification_homer_library_type_gene_table:
    input:
        gq_homer_prefix + 'index.{norm_method}.{library_type}.txt'
    output:
        gq_homer_prefix + 'gene.{norm_method}.{library_type}.tsv'
    run:
        combined_df = combine_quantification_tables(index_filepath = input[0], 
                                                    output_filepath = output[0], 
                                                    df_filepath_suffix = f'{wildcards.norm_method}.tsv', 
                                                    key_column = 'gene_id', 
                                                    value_column = wildcards.norm_method, 
                                                    renamed_key_column = 'gene_id',
                                                    sep = '\t')


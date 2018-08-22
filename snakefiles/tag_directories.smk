from itertools import compress

# Import tag directory default configuration
configfile: 'configs/tag_directories.yaml'

# Import tag directory configuration as variables
tag_directories_prefix = config['tag_directories']['prefix']
tag_directories_library_types = config['tag_directories']['library_types']
# We assume only paired end and single end libraries for now
end_types = ['paired_end', 'single_end']
tag_directories_input_prefix_by_library_type = {}
tag_directories_options_by_library_type = {}
for library_type in tag_directories_library_types:
    tag_directories_input_prefix_by_library_type[library_type] = {end_type: tag_directories_library_types[library_type][end_type]['input_prefix'] for end_type in end_types}
    tag_directories_options_by_library_type[library_type] = {end_type: tag_directories_library_types[library_type][end_type]['params']['options'] for end_type in end_types}

# Attempt paired-end options first, but if a file containing the second mate reads isn't available, then go for single end options
ruleorder: tag_directories_homer_paired_end > tag_directories_homer_single_end

rule tag_directories_homer_paired_end:
    input:
        bam = lambda wildcards: (tag_directories_input_prefix_by_library_type[wildcards.library_type]['paired_end'] + basename_pattern + '.sorted.bam').format(**wildcards),
        reads2 = reads_prefix + basename_pattern + reads2_suffix + fastq_suffix + compression_suffix
    output:
        tag_directories_prefix + basename_pattern + '/tagInfo.txt'
    log:
        tag_directories_prefix + basename_pattern + '.log'
    params:
        options = lambda wildcards: tag_directories_options_by_library_type[wildcards.library_type]['paired_end']
    shadow:
        'shallow'
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        makeTagDirectory $(dirname {output}) {params.options} {input.bam} &> {log}
        """

rule tag_directories_homer_single_end:
    input:
        bam = lambda wildcards: (tag_directories_input_prefix_by_library_type[wildcards.library_type]['single_end'] + basename_pattern + '.sorted.bam').format(**wildcards)
    output:
        tag_directories_prefix + basename_pattern + '/tagInfo.txt'
    log:
        tag_directories_prefix + basename_pattern + '.log'
    params:
        options = lambda wildcards: tag_directories_options_by_library_type[wildcards.library_type]['single_end']
    shadow:
        'shallow'
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        makeTagDirectory $(dirname {output}) {params.options} {input.bam} &> {log}
        """

rule tag_directories_homer_index_from_reads:
    input:
        lambda wildcards: [(tag_directories_prefix + basename + '/tagInfo.txt').format(**wildcards) for basename in get_read_basenames_for_organism(wildcards.organism)]
    output:
        tag_directories_prefix + 'index.txt'
    run:
        shell(f'mkdir -p $(dirname {output[0]})')
        tag_info_list = [tag_info for tag_info in input]
        tag_dirs_list = sorted(['/'.join(tag_info.split('/')[:-1])+'/' for tag_info in tag_info_list])
        tag_dirs_string = '\n'.join(tag_dirs_list)+'\n'
        with open(output[0], 'w') as outfile:
            outfile.write(tag_dirs_string)

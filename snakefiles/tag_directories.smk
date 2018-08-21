# Import naming convention for reads
configfile: 'configs/reads.yaml'
# Import tag directory default configuration
configfile: 'configs/tag_directories.yaml'

# Import naming conventions for reads as variables
reads_prefix = config['reads']['prefix']
reads1_suffix = config['reads']['reads1_suffix']
reads2_suffix = config['reads']['reads2_suffix']
fastq_suffix = config['reads']['fastq_suffix']
compression_suffix = config['reads']['compression_suffix']

# Import tag directory configuration as variables
## Basename pattern contains wildcards in curly brace format that are referenced. the most important one is the {library_type} wildcard
basename_pattern = config['tag_directories']['basename_pattern']
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

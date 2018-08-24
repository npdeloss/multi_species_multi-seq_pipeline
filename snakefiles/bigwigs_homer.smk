# Import naming conventions for reads as variables
reads_prefix = config['reads']['prefix']
reads1_suffix = config['reads']['reads1_suffix']
reads2_suffix = config['reads']['reads2_suffix']
fastq_suffix = config['reads']['fastq_suffix']
compression_suffix = config['reads']['compression_suffix']

# Import naming conventions for bigwig generation as variables
bigwigs_homer_prefix = config['bigwigs']['homer']['prefix']
bigwigs_homer_input_prefix = config['bigwigs']['homer']['input_prefix']
bigwigs_homer_genome_prefix = config['bigwigs']['homer']['genome_prefix']
bigwigs_homer_options = config['bigwigs']['homer']['params']['options']

ruleorder: bigwigs_homer > bigwigs_homer_alias_both

rule bigwigs_homer:
    input:
        tag_info = bigwigs_homer_input_prefix + '{basename}/tagInfo.txt',
        chrom_sizes = bigwigs_homer_genome_prefix + 'genome.chrom.sizes'
    output:
        bigwig = bigwigs_homer_prefix + '{basename}.{strand}.bw',
        track_info = bigwigs_homer_prefix + '{basename}.{strand}.trackInfo.txt'
    log:
        bigwig = bigwigs_homer_prefix + '{basename}.{strand}.log'
    params:
        options = bigwigs_homer_options,
        strand = lambda wildcards: {'both':'both','pos':'+', 'neg':'-'}[wildcards.strand]
    shadow:
        'shallow'
    conda:
        '../envs/homer.yaml'
    shell:
        """
        mkdir -p $(dirname {output.bigwig})
        makeUCSCfile \
        $(dirname {input.tag_info}) \
        -o {output.bigwig} \
        -bigWig {input.chrom_sizes} \
        {params.options} \
        -strand {params.strand} \
        > {output.track_info} \
        2> {log}
        """

rule bigwigs_homer_alias_both:
    input:
        bigwig = bigwigs_homer_prefix + '{basename}.both.bw',
        track_info = bigwigs_homer_prefix + '{basename}.both.trackInfo.txt'
    output:
        bigwig = bigwigs_homer_prefix + '{basename}.bw',
        track_info = bigwigs_homer_prefix + '{basename}.trackInfo.txt'
    shell:
        """
        mkdir -p $(dirname {output.bigwig})
        ln -sf $(basename {input.bigwig}) {output.bigwig}
        ln -sf $(basename {input.track_info}) {output.track_info}
        """

rule bigwigs_homer_index_library_type:
    input:
        lambda wildcards: [(bigwigs_homer_prefix + basename + '.{strand}.bw').format(**wildcards) for basename in get_read_basenames_for_organism(wildcards.organism, 
                                                                                                                                                        library_type = wildcards.library_type)]
    output:
        bigwigs_homer_prefix + 'index.{strand}.{library_type}.txt'
    run:
        shell(f'mkdir -p $(dirname {output[0]})')
        input_files = [input_file for input_file in input]
        input_files_string = '\n'.join(input_files)+'\n'
        with open(output[0], 'w') as outfile:
            outfile.write(input_files_string)

rule bigwigs_homer_index_library_type_alias_both:
    input:
        lambda wildcards: [(bigwigs_homer_prefix + basename + '.bw').format(**wildcards) for basename in get_read_basenames_for_organism(wildcards.organism, 
                                                                                                                                                        library_type = wildcards.library_type)]
    output:
        bigwigs_homer_prefix + 'index.{library_type}.txt'
    run:
        shell(f'mkdir -p $(dirname {output[0]})')
        input_files = [input_file for input_file in input]
        input_files_string = '\n'.join(input_files)+'\n'
        with open(output[0], 'w') as outfile:
            outfile.write(input_files_string)

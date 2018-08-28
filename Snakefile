import os

configfile: 'configs/wildcard_constraints.yaml'
wildcard_constraints:
    **config['wildcard_constraints']

# Import naming convention for reads
configfile: 'configs/reads.yaml'

# Import naming conventions for reads as variables
## Basename pattern contains wildcards in curly brace format that are referenced. The most important ones are organism and library_type wildcards.
basename_pattern = config['reads']['basename_pattern']
reads_prefix = config['reads']['prefix']
reads1_suffix = config['reads']['reads1_suffix']
reads2_suffix = config['reads']['reads2_suffix']
fastq_suffix = config['reads']['fastq_suffix']
compression_suffix = config['reads']['compression_suffix']

# Get all basenames for reads in our reads directory that match a particular organism (and optionally, library type)
def get_read_basenames_for_organism(organism, library_type = None, 
                                    basename_pattern = basename_pattern, 
                                    reads_prefix = reads_prefix, 
                                    suffix = reads1_suffix + fastq_suffix + compression_suffix):
    # Concatenate the pattern to search for reads
    reads_pattern = reads_prefix + basename_pattern + suffix
    # Search for files, get wildcard values
    wc_dict = dict(glob_wildcards(reads_pattern)._asdict())
    # If we aren't looking for any organism in particular, just return all the basenames
    if organism is None:
        return sorted(expand(basename_pattern, zip, **wc_dict))
    # Filter only for wildcard values at indices where the organism matches our argument
    organism_index = [organism_entry == organism for organism_entry in wc_dict['organism']]
    wc_dict_for_organism = {key: list(compress(wc_dict[key], organism_index)) for key in wc_dict}
    # If we aren't asking for a specific libary type, expand the basenames out with the filtered wildcards and return them
    if library_type is None:
        return sorted(expand(basename_pattern, zip, **wc_dict_for_organism))
    # Otherwise, filter for wildcards with matching library types too
    library_type_index = [library_type_entry == library_type for library_type_entry in wc_dict_for_organism['library_type']]
    wc_dict_for_library_type = {key: list(compress(wc_dict_for_organism[key], library_type_index)) for key in wc_dict_for_organism}
    # Expand the basenames out with the filtered wildcards and return them
    return sorted(expand(basename_pattern, zip, **wc_dict_for_library_type))

include: 'snakefiles/genomes.smk'
include: 'snakefiles/preprocessing.smk'
include: 'snakefiles/qc.smk'
include: 'snakefiles/alignment.smk'
include: 'snakefiles/tag_directories.smk'
include: 'snakefiles/bigwigs.smk'
include: 'snakefiles/gene_quantifications.smk'
include: 'snakefiles/track_visualizations.smk'

configfile: 'config.yaml'

rule targets_local:
    input:
        *config['targets']
    output:
        'targets_local.txt'
    run:
        with open(output[0],'w') as outfile:
            outfile.write('\n'.join([inputfile for inputfile in input])+'\n')

rule www:
    input:
        targets_local = 'targets_local.txt', 
        dirs = lambda wildcards: [directory(config['www']['dir'] + source_dir) for source_dir in config['www']['source_dirs']],
        html = 'igv_js.html'
    output:
        igv_js = config['www']['dir'] + 'igv_js.html',
        www = 'www'
    shell:
        """
        mkdir -p $(dirname {output.igv_js})
        cp {input.html} {output.igv_js}
        cp *.html $(dirname {output.igv_js})
        touch {output.www}
        """

rule www_from_source_dir:
    input:
        *config['targets']
    output:
        directory(config['www']['dir'] + '{source_dir}')
    log:
        '{source_dir}.www_rsync.log'
    run:
        if os.path.isdir(wildcards.source_dir):
            shell(f'mkdir -p {config["www"]["dir"]}')
            shell(f'rsync  --verbose --archive --recursive {wildcards.source_dir} {config["www"]["dir"]} &> {log[0]}')

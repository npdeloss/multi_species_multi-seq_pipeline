configfile: 'configs/wildcard_constraints.yaml'
wildcard_constraints:
    **config['wildcard_constraints']

include: 'snakefiles/genomes.smk'
include: 'snakefiles/preprocessing.smk'
include: 'snakefiles/qc.smk'
include: 'snakefiles/alignment.smk'
include: 'snakefiles/tag_directories.smk'
include: 'snakefiles/bigwigs.smk'

configfile: 'config.yaml'

rule all:
    input:
        **config['targets']

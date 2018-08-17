include: 'snakefiles/genomes.smk'
include: 'snakefiles/preprocessing.smk'
include: 'snakefiles/qc.smk'
include: 'snakefiles/alignment.smk'

configfile: 'config.yaml'

rule all:
    input:
        **config['targets']

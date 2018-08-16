include: 'genomes.rules'

configfile: 'config.yaml'

rule all:
    input:
        **config['targets']

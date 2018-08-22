# Import quantification default configuration
configfile: 'configs/gene_quantifications.yaml'

# Include different quantification methods
include: 'gene_quantifications_kallisto.smk'
# include: 'gene_quantifications_homer.smk'

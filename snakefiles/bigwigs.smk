# Import naming convention for reads
configfile: 'configs/reads.yaml'
# Import bigwig default configuration
configfile: 'configs/bigwigs.yaml'

# Include different bigwig generators
include: 'bigwigs_homer.smk'

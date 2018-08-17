# Import naming convention for reads
configfile: 'configs/reads.yaml'
# Import alignment default configuration
configfile: 'configs/alignment.yaml'

include: 'alignment_star.smk'
# include: 'alignment_bowtie2.smk'

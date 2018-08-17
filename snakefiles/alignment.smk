# Import naming convention for reads
configfile: 'configs/reads.yaml'
# Import alignment default configuration
configfile: 'configs/alignment.yaml'

include: 'snakefiles/alignment_star.smk'
include: 'snakefiles/alignment_bowtie2.smk'
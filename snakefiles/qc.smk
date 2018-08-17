# Import naming conventions for reads
configfile: 'configs/reads.yaml'
# Import preprocessing QC configuration
configfile: 'configs/qc.yaml'

# Import user configuration
configfile: 'config.yaml'

include: 'snakefiles/qc_fastqc.smk'
# Import genome default configuration
configfile: 'configs/genomes.yaml'

# Import user configuration
configfile: 'config.yaml'

prefix = config['genomes']['prefix']

# Download and decompress reference files
rule download_genome:
    output:
        prefix+'{organism}/{reference}/downloads/genome.fa.gz'
    params:
        url = lambda wildcards: config['genomes'][wildcards.organism][wildcards.reference]['genome_fa']['url']
    log:
        prefix+'{organism}/{reference}/downloads/genome.fa.gz.log'
    conda:
        '../envs/wget.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        wget -O {output} -o {log} {params.url}
        """

rule download_annotation:
    output:
        prefix+'{organism}/{reference}/downloads/annotation.gtf.gz'
    params:
        url = lambda wildcards: config['genomes'][wildcards.organism][wildcards.reference]['annotation_gtf']['url']
    log:
        'genomes/{organism}/{reference}/downloads/annotation.gtf.gz.log'
    conda:
        '../envs/wget.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        wget -O {output} -o {log} {params.url}
        """

rule decompress_reference_file:
    input:
        prefix+'{organism}/{reference}/downloads/{basename}.gz'
    output:
        prefix+'{organism}/{reference}/{basename}'
    shell:
        """
        zcat {input} > {output}
        """

# Index Genome
rule index_genome:
    input:
        prefix+'{organism}/{reference}/genome.fa'
    output:
        prefix+'{organism}/{reference}/genome.fa.fai'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        samtools faidx {input}
        """

rule chrom_size_genome:
    input:
        prefix+'{organism}/{reference}/genome.fa.fai'
    output:
        prefix+'{organism}/{reference}/genome.chrom.sizes'
    shell:
        """
        cut -f1,2 {input} > {output}
        """

# Index Annotation
rule bgzip_annotation:
    input:
        prefix+'{organism}/{reference}/annotation.gtf'
    output:
        prefix+'{organism}/{reference}/annotation.gtf.gz'
    log:
        prefix+'{organism}/{reference}/annotation.gtf.gz.log'
    conda:
        '../envs/tabix.yaml'
    shell:
        """
        (cat {input} \
        | grep -v '^#' \
        | sort -k1,1 -k4,4n \
        | bgzip -c \
        > {output}) \
        2> {log}
        """

rule tabix_annotation:
    input:
        prefix+'{organism}/{reference}/annotation.gtf.gz'
    output:
        prefix+'{organism}/{reference}/annotation.gtf.gz.tbi'
    log:
        prefix+'{organism}/{reference}/annotation.gtf.gz.tbi.log'
    conda:
        '../envs/tabix.yaml'
    shell:
        """
        sleep 5s
        tabix -f -p gff {input} 2> {log}
        """

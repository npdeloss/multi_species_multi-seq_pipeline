configfile: 'configs/wildcard_constraints.yaml'
wildcard_constraints:
    **config['wildcard_constraints']

# Import genome default configuration
configfile: 'configs/genomes.yaml'

# Import user configuration
configfile: 'config.yaml'

prefix = config['genomes']['prefix']

# Download and decompress reference files
rule download_genome:
    output:
        prefix+'downloads/genome.fa.gz'
    params:
        url = lambda wildcards: config['genomes'][wildcards.organism][wildcards.reference]['genome_fa']['url']
    log:
        prefix+'downloads/genome.fa.gz.log'
    conda:
        '../envs/wget.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        wget -O {output} -o {log} {params.url}
        """

rule download_annotation:
    output:
        prefix+'downloads/annotation.gtf.gz'
    params:
        url = lambda wildcards: config['genomes'][wildcards.organism][wildcards.reference]['annotation_gtf']['url']
    log:
        prefix+'downloads/annotation.gtf.gz.log'
    conda:
        '../envs/wget.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        wget -O {output} -o {log} {params.url}
        """

rule decompress_genoma_fa:
    input:
        prefix+'downloads/genome.fa.gz'
    output:
        prefix+'genome.fa'
    shell:
        """
        zcat {input} > {output}
        """

rule decompress_annotation_gtf:
    input:
        prefix+'downloads/annotation.gtf.gz'
    output:
        prefix+'annotation.gtf'
    shell:
        """
        zcat {input} > {output}
        """

# Index Genome
rule index_genome:
    input:
        prefix+'genome.fa'
    output:
        prefix+'genome.fa.fai'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        samtools faidx {input}
        """

rule chrom_size_genome:
    input:
        prefix+'genome.fa.fai'
    output:
        prefix+'genome.chrom.sizes'
    shell:
        """
        cut -f1,2 {input} > {output}
        """

# Index Annotation
rule bgzip_annotation:
    input:
        prefix+'annotation.gtf'
    output:
        prefix+'annotation.gtf.gz'
    log:
        prefix+'annotation.gtf.gz.log'
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
        prefix+'annotation.gtf.gz'
    output:
        prefix+'annotation.gtf.gz.tbi'
    log:
        prefix+'annotation.gtf.gz.tbi.log'
    conda:
        '../envs/tabix.yaml'
    shell:
        """
        sleep 5s
        tabix -f -p gff {input} 2> {log}
        """

rule annotation_tsv:
    input:
        prefix+'annotation.gtf'
    output:
        prefix+'annotation.tsv'
    log:
        prefix+'annotation.tsv.log'
    conda:
        '../envs/gtfparse.yaml'
    shell:
        """
        python scripts/gtf_to_tsv.py {input} {output} &> {log}
        """

rule annotation_simplified_bed:
    input:
        prefix+'annotation.tsv'
    output:
        prefix+'annotation.simplified.bed'
    run:
        annotation = pd.read_table(input[0]).query('feature == "gene"')
        annotation_bed = annotation[['seqname', 'start', 'end', 'gene_name', 'score', 'strand']].drop_duplicates().sort_values(['seqname', 'start'])
        annotation_bed['score'] = 0
        annotation_bed.to_csv(output[0], index = False, header = False)

rule generate_transcriptome:
    input:
        genome_fa = prefix+'genome.fa',
        annotation_gtf = prefix+'annotation.gtf'
    output:
        prefix+'transcriptome.fa'
    log:
        prefix+'transcriptome.log'
    conda:
        '../envs/gffread.yaml'
    shell:
        """
        gffread -w {output} -g {input.genome_fa} {input.annotation_gtf} &> {log}
        """

rule index_transcriptome:
    input:
        prefix+'transcriptome.fa'
    output:
        prefix+'transcriptome.fa.fai'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        samtools faidx {input}
        """

rule seq_size_transcriptome:
    input:
        prefix+'transcriptome.fa.fai'
    output:
        prefix+'transcriptome.seq.sizes'
    shell:
        """
        cut -f1,2 {input} > {output}
        """

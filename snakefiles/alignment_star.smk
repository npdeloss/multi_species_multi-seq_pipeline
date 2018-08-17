# Import naming conventions for reads as variables
reads_prefix = config['reads']['prefix']
reads1_suffix = config['reads']['reads1_suffix']
reads2_suffix = config['reads']['reads2_suffix']
fastq_suffix = config['reads']['fastq_suffix']
compression_suffix = config['reads']['compression_suffix']

# Import prefixes for input and output
## Indexing
index_star_prefix = config['alignment']['star']['index']['prefix']
index_star_input_prefix = config['alignment']['star']['index']['input_prefix']
index_star_threads = config['alignment']['star']['index']['threads']
index_star_params_options = config['alignment']['star']['index']['params']['options']

## Alignment
alignment_star_threads = config['alignment']['star']['threads']
alignment_star_prefix = config['alignment']['star']['prefix']
alignment_star_input_prefix = config['alignment']['star']['input_prefix']
alignment_star_read_files_command = config['alignment']['star']['read_files_command']

### Paired end
alignment_star_paired_end_params_options = config['alignment']['star']['paired_end']['params']['options']
alignment_star_single_end_params_options = config['alignment']['star']['single_end']['params']['options']

# Index genome
rule index_star:
    input:
        genome_fa = index_star_input_prefix + 'genome.fa',
        annotation_gtf = index_star_input_prefix + 'annotation.gtf',
    output:
        alignment_star_prefix + 'Genome'
    log:
        alignment_star_prefix + 'index_star.log'
    threads:
        index_star_threads
    params:
        options = index_star_params_options
    conda:
        '../envs/star.yaml'
    shell:
        """
        mkdir -p  $(dirname {output})
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome_fa} \
        --sjdbGTFfile {input.annotation_gtf} \
        --outFileNamePrefix {output}/ \
        --outStd Log &> {log}
        """

# Attempt paired end alignment before single end alignment
ruleorder: alignment_star_paired_end > alignment_star_single_end

# Paired end alignment
rule alignment_star_paired_end:
    input:
        reads1 = alignment_star_input_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        reads2 = alignment_star_input_prefix + '{basename}' + reads2_suffix + fastq_suffix + compression_suffix,
        index = alignment_star_prefix + 'Genome'
    output:
        bam = alignment_star_prefix + '{basename}/Aligned.out.bam',
        counts = alignment_star_prefix + '{basename}/ReadsPerGene.out.tab'
    log:
        alignment_star_prefix + '{basename}.log'
    threads:
        alignment_star_threads
    params:
        read_files_command = alignment_star_read_files_command,
        options = alignment_star_paired_end_params_options
    conda:
        '../envs/star.yaml'
    shell:
        """
        outdir=$(dirname {output.bam})
        mkdir -p $outdir
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.index} \
        --readFilesIn {input.reads1} {input.reads2} \
        --readFilesCommand {params.read_files_command} \
        --outSAMtype BAM Unsorted \
        --quantMode GeneCounts \
        --outFileNamePrefix $outdir/ \
        {params.options} \
        --outStd Log &> {log}
        """

# Single end alignment
rule alignment_star_single_end:
    input:
        reads1 = alignment_star_input_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        index = alignment_star_prefix + 'Genome'
    output:
        bam = alignment_star_prefix + '{basename}/Aligned.out.bam',
        counts = alignment_star_prefix + '{basename}/ReadsPerGene.out.tab'
    log:
        alignment_star_prefix + '{basename}.log'
    threads:
        alignment_star_threads
    params:
        read_files_command = alignment_star_read_files_command,
        options = alignment_star_single_end_params_options
    conda:
        '../envs/star.yaml'
    shell:
        """
        outdir=$(dirname {output.bam})
        mkdir -p $outdir
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.index} \
        --readFilesIn {input.reads1} \
        --readFilesCommand {params.read_files_command} \
        --outSAMtype BAM Unsorted \
        --quantMode GeneCounts \
        --outFileNamePrefix $outdir/ \
        {params.options} \
        --outStd Log &> {log}
        """

# Sort alignment
alignment_sort_threads = config['alignment']['sort']['threads']
alignment_sort_params_options = config['alignment']['sort']['params']['options']

rule alignment_star_sort_bam:
    input:
        alignment_star_prefix + '{basename}/Aligned.out.bam'
    output:
        alignment_star_prefix + '{basename}.sorted.bam'
    threads:
        alignment_sort_threads
    log:
        alignment_star_prefix + '{basename}.sorted.log'
    params:
        options = alignment_sort_params_options
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        samtools sort \
        -o {output} \
        -T {output}.sorting.tmp \
        -@ {threads} \
        {params.options} \
        {input} &> {log}
        """

# Index sorted alignment
rule alignment_star_index_bam:
    input:
        alignment_star_prefix + '{basename}.sorted.bam'
    output:
        alignment_star_prefix + '{basename}.sorted.bam.bai'
    log:
        alignment_star_prefix + '{basename}.sorted.index.log'
    shell:
        """
        sleep 1s
        samtools index {input} &> {log}
        """

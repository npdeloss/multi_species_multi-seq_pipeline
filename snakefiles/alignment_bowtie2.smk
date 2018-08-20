# Import naming conventions for reads as variables
reads_prefix = config['reads']['prefix']
reads1_suffix = config['reads']['reads1_suffix']
reads2_suffix = config['reads']['reads2_suffix']
fastq_suffix = config['reads']['fastq_suffix']
compression_suffix = config['reads']['compression_suffix']

# Import configuration options
## Indexing
index_bowtie2_prefix = config['alignment']['bowtie2']['index']['prefix']
index_bowtie2_input_prefix = config['alignment']['bowtie2']['index']['input_prefix']
index_bowtie2_threads = config['alignment']['bowtie2']['index']['threads']
index_bowtie2_params_options = config['alignment']['bowtie2']['index']['params']['options']

## Alignment
alignment_bowtie2_threads = config['alignment']['bowtie2']['threads']
alignment_bowtie2_prefix = config['alignment']['bowtie2']['prefix']
alignment_bowtie2_input_prefix = config['alignment']['bowtie2']['input_prefix']

### Paired end
alignment_bowtie2_paired_end_params_options = config['alignment']['bowtie2']['paired_end']['params']['options']

## Single end
alignment_bowtie2_single_end_params_options = config['alignment']['bowtie2']['single_end']['params']['options']

## Sorting alignment
alignment_sort_threads = config['alignment']['sort']['threads']
alignment_sort_params_options = config['alignment']['sort']['params']['options']

# Index genome
rule index_bowtie2:
    input:
        genome_fa = index_bowtie2_input_prefix + 'genome.fa'
    output:
        index_bowtie2_prefix + 'index.1.bt2'
    log:
        index_bowtie2_prefix + 'index_bowtie2.log'
    threads:
        index_bowtie2_threads
    params:
        options = index_bowtie2_params_options
    conda:
        '../envs/bowtie2.yaml'
    shell:
        """
        mkdir -p $(dirname {output})
        bowtie2-build \
        {params.options} \
        {input.genome_fa} \
        $(dirname {output})/$(basename {output} .1.bt2) \
        &> {log}
        """

# Attempt paired end alignment before single end alignment
ruleorder: alignment_bowtie2_paired_end > alignment_bowtie2_single_end

# Paired end alignment
rule alignment_bowtie2_paired_end:
    input:
        reads1 = alignment_bowtie2_input_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        reads2 = alignment_bowtie2_input_prefix + '{basename}' + reads2_suffix + fastq_suffix + compression_suffix,
        index = index_bowtie2_prefix + 'index.1.bt2'
    output:
        bam = alignment_bowtie2_prefix + '{basename}/Aligned.out.bam'
    log:
        alignment_bowtie2_prefix + '{basename}.log'
    threads:
        alignment_bowtie2_threads
    params:
        options = alignment_bowtie2_paired_end_params_options
    conda:
        '../envs/bowtie2.yaml'
    shell:
        """
        out_dir=$(dirname {output.bam})
		index=$(dirname {input.index})/$(basename {input.index} .1.bt2)
        mkdir -p $out_dir
        (bowtie2 \
		{params.options} \
        -p {threads} \
        -x $index \
        -1 {input.reads1} \
        -2 {input.reads2} \
        | samtools sort \
        -n \
        -O bam \
        -T $out_dir/.sorting \
        -@ {threads} - \
        | samtools fixmate -m - - \
        | samtools sort - -@ {threads} \
        | samtools markdup - {output.bam}) 2> {log}
        """

# Single end alignment
rule alignment_bowtie2_single_end:
    input:
        reads1 = alignment_bowtie2_input_prefix + '{basename}' + reads1_suffix + fastq_suffix + compression_suffix,
        index = index_bowtie2_prefix + 'index.1.bt2'
    output:
        bam = alignment_bowtie2_prefix + '{basename}/Aligned.out.bam'
    log:
        alignment_bowtie2_prefix + '{basename}.log'
    threads:
        alignment_bowtie2_threads
    params:
        options = alignment_bowtie2_paired_end_params_options
    conda:
        '../envs/bowtie2.yaml'
    shell:
        """
        out_dir=$(dirname {output.bam})
		index=$(dirname {input.index})/$(basename {input.index} .1.bt2)
        mkdir -p $out_dir
        (bowtie2 \
		{params.options} \
        -p {threads} \
        -x $index \
        -U {input.reads1} \
        | samtools sort \
        -n \
        -O bam \
        -T $out_dir/.sorting \
        -@ {threads} - \
        | samtools fixmate -m - - \
        | samtools sort - -@ {threads} \
        | samtools markdup - {output.bam}) 2> {log}
		"""

# Sort alignment

rule alignment_bowtie2_sort_bam:
    input:
        alignment_bowtie2_prefix + '{basename}/Aligned.out.bam'
    output:
        alignment_bowtie2_prefix + '{basename}.sorted.bam'
    threads:
        alignment_sort_threads
    log:
        alignment_bowtie2_prefix + '{basename}.sorted.log'
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
rule alignment_bowtie2_index_bam:
    input:
        alignment_bowtie2_prefix + '{basename}.sorted.bam'
    output:
        alignment_bowtie2_prefix + '{basename}.sorted.bam.bai'
    log:
        alignment_bowtie2_prefix + '{basename}.sorted.index.log'
    conda:
        '../envs/samtools.yaml'
    shell:
        """
        sleep 1s
        samtools index {input} &> {log}
        """-

# Snakemake and Conda -based RNA-seq and ChIP-seq Pipeline
This is a snakefile and accompanying conda environments for processing the Benner Lab's ChIP-seq and RNA-seq data. It adapts work by Carlos Guzman and Samuel Roth into a directory structure that supports multiple species and reference annotations in a single project. It relies on a uniform input file naming scheme to automate tasks such as peak calling and alignment after deposition of all project fast files into a single folder.

## Getting Started
These instructions will get you a copy of the project up and running on your local machine. As a goal, all configuration should be possible without privilege escalation.

### Prerequisites
To run this snakefile, you need only Snakemake and Conda/Miniconda. Future releases may require pandas for parsing datasheets.
While the pipeline itself relies on many software packages, most notable Homer, these dependencies are handled using Snakemake's ability to specify and download conda environments for individual rules. This allows rules to depend on their own versions of packages, rather than relying on all package versions being compatible with each other across all rules.

Instructions on installing Conda can be found on its [Conda documentation page](https://conda.io/docs/user-guide/install/index.html).  
After installing Conda, Snakemake can be installed with Conda as follows:
```
conda install -c bioconda snakemake 
```

### Installing
To install, copying `Snakefile`, `config.yaml`, and the `envs` folder from this repository will suffice. Put them in the same folder.

### Usage
#### Naming and depositing fastq files
First, ensure your fastq files follow the appropriate naming scheme. The naming scheme should be as follows:  
`[human/mouse]_{cell}_[rna/chip]_{description}_{initials}_[R1/R2]_001.fastq.gz`  
An example is:  
`human_htbe_chip_h3k27ac-cond1-rep1_YC_R1_001.fastq.gz`  
A few notes on these fields:
* Currently, only human and mouse genomes are supported via GENCODE annotation.
* Values for `{cell}` and `{libary_subtype}` are restricted to only be alphanumeric.
* The `{description}` may not include spaces or underscores.
* The `{initials}` field refers to the individual who performed the experiment. This should be alphanumeric only.
* Paired-end sequencing is detected by finding files with both `_R1_001.fastq.gz` and `_R2_001.fastq.gz` suffixes.

If your data are not named in this manner, you can create symlinks to them like so:  
`ln -sf [old_filename] [new_filename]`  

To deposit files into this project for use with the pipeline, create a folder called `reads_fastq` in the same directory as the Snakefile and copy them into there (symlinking is also an option).  
```
mkdir -p fastq_files
cp /your/directory/here/*.fastq.gz reads_fastq/
```

### Calling Snakemake
To use snakemake's dependency downloading capability through conda, invoke the pipeline using:
```
snakemake --use-conda -j {threads} {target}
```
Where `threads` is the number of threads you wish to use, and `{target}` is the name of the target file. More on these later.

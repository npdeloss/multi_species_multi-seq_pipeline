# Snakemake and Conda -based RNA-seq and ChIP-seq Pipeline
This is a snakefile and accompanying conda environments for processing the Benner Lab's ChIP-seq and RNA-seq data. It adapts work by Carlos Guzman and Samuel Roth into a directory structure that supports multiple species and reference annotations in a single project. It relies on a uniform input file naming scheme to automate tasks such as peak calling and alignment after deposition of all project fast files into a single folder.

## Getting Started
These instructions will get you a copy of the project up and running on your local machine. As a goal, all configuration should be possible without privilege escalation.

### Prerequisites
To run this snakefile, you need only Snakemake and Conda/Miniconda. Future releases may require pandas for parsing datasheets.
While the pipeline itself relies on many software packages, most notably [HOMER](http://homer.ucsd.edu/homer/), these dependencies are handled using Snakemake's ability to specify and download conda environments for individual rules. This allows rules to depend on their own versions of packages, rather than relying on all package versions being compatible with each other across all rules.

Instructions on installing Conda can be found on its [Conda documentation page](https://conda.io/docs/user-guide/install/).  
After installing Conda, Snakemake can be installed with Conda as follows:
```
conda install -c bioconda snakemake 
```

### Installing
To install, copying `Snakefile`, `config.yaml`, `igv_js.html`, `envs`, `configs`, `resources`, `scripts`, and `snakefiles` from this repository will suffice. Put them in the same folder. Alternatively, download this repository and rename as required.

```
git clone https://github.com/npdeloss/multi_species_multi-seq_pipeline.git
mv multi_species_multi-seq_pipeline your_project_name
cd your_project_name
```

### Usage
#### Naming and depositing fastq files
First, ensure your fastq files follow the appropriate naming scheme. The naming scheme should be as follows:  
`[human/mouse]_{cell}_[rna/chip]_{description}_{initials}_[R1/R2]_001.fastq.gz`  
An example is:  
`human_htbe_chip_h3k27ac-cond1-rep1_YC_R1_001.fastq.gz`  
A few notes on these fields:
* Currently, only human and mouse genomes are supported via GENCODE annotation.
* Values for `{cell}` and `{libary_type}` may not contain underscores or slashes.
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

#### Calling Snakemake
To use snakemake's dependency downloading capability through conda, invoke the pipeline using:
```
snakemake --use-conda -j {threads} {target}
```
Where `threads` is the number of threads you wish to use, and `{target}` is the name of the target file. More on these later.

#### Quick targets
To run all available analyses automatically, run:
```
snakemake --use-conda -j {threads} targets
```

You may need to remove the `targets` file if it has already been generated:
```
rm targets
```

#### Web browser targets
Change the second line of `config.yaml` to reflect the directory where your server stores files on your project to serve to web browsers. This is usually `/var/www/html/your_project_name`. Change the third line to reflect the URL where your project will be served to the web. Usually, this is `http://your.domain.com/your_project_name`.

```
www:
    dir: '/var/www/html/your_project_name/'
    url: 'http://your.domain.com/your_project_name/'
```

#### Targets in detail
In `config.yaml`, the `targets` object contains lists of targets that will be created by the pipeline by default when running analyses automatically. These are decribed within config.yaml for the `human/grch38` organism/reference combination.

##### Generate FASTQC reports for all samples
If the library_type is chip, align with bowtie2. If library_type is rna, align with STAR.
```
qc_fastqc/index.txt
```

##### Gene quantifications using kallisto
Quantifies gene expression for all samples with `rna` `library_type`.
```
# TPM for basic normalization
gene_quantifications_kallisto/human/grch38/gene.tpm.rna.tsv
# Estimated counts, for use with programs expecting count-like data.
gene_quantifications_kallisto/human/grch38/gene.est_counts.rna.tsv
```
##### Gene quantifications using HOMER for all samples with rna library_type
Quantifies gene expression for all samples with `rna` `library_type`.
```
# Read counts
gene_quantifications_homer/human/grch38/gene.raw.rna.tsv
# Reads normalized to 1e7 by default
gene_quantifications_homer/human/grch38/gene.norm.rna.tsv
# Expression normalized with RPKM
gene_quantifications_homer/human/grch38/gene.rpkm.rna.tsv
```
##### Bigwig coverage tracks using HOMER
For all samples with rna library_type
```
# Positive strand
bigwigs_homer/human/grch38/index.pos.rna.txt
# Negative strand
bigwigs_homer/human/grch38/index.neg.rna.txt
# Both strands
bigwigs_homer/human/grch38/index.rna.txt
```

For all samples with chip library_type
```
bigwigs_homer/human/grch38/index.chip.txt
```

##### JSON and HTML for visualizing the bigwig coverage tracks using IGV.js
These only work when served from a server. Paths are relative to the `dir` property of the `www` section from `config.yaml`
For example, if the `dir` property was set to `/var/www/html/your_project_name/`, and the `url` propery was set to `http://your.domain.com/your_project_name/`, then you would point your browser to `http://your.domain.com/your_project_name/track_visualizations_bigwigs/human/grch38/chip.igv.html` to visualize coverage for your ChIP-seq tracks.
```
# JSON track data
track_visualizations_bigwigs/human/grch38/rna.igv.json
track_visualizations_bigwigs/human/grch38/chip.igv.json
# HTML page
track_visualizations_bigwigs/human/grch38/rna.igv.html
track_visualizations_bigwigs/human/grch38/chip.igv.html
```

# Import quantification default configuration
configfile: 'configs/gene_quantifications.yaml'
import pandas as pd

def combine_quantification_tables(index_filepath, output_filepath, df_filepath_suffix, key_column, value_column, sep = '\t'):
    dirs = [line.strip() for line in open(index_filepath).readlines()]
    dfs_by_sample = {dirname.split('/')[-2]: pd.read_table(dirname+df_filepath_suffix).set_index(key_column) for dirname in dirs}
    samples = sorted(list(dfs_by_sample.keys()))
    series_by_sample = {sample: dfs_by_sample[sample][value_column].rename(sample) for sample in samples}
    combined_df = pd.DataFrame(series_by_sample).reset_index()
    combined_df.to_csv(output_filepath, sep = sep, index = False)
    return combined_df

# Include different quantification methods
include: 'gene_quantifications_kallisto.smk'
# include: 'gene_quantifications_homer.smk'

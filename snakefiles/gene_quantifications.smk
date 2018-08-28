# Import quantification default configuration
configfile: 'configs/gene_quantifications.yaml'
import pandas as pd

def combine_quantification_tables(index_filepath, output_filepath, df_filepath_suffix, key_column, value_column, renamed_key_column = None, sep = '\t'):
    dirs = [line.strip() for line in open(index_filepath).readlines() if line.strip() is not '']
    if len(dirs) > 0:
        dfs_by_sample = {dirname.split('/')[-2]: pd.read_table(dirname+df_filepath_suffix).set_index(key_column) for dirname in dirs}
        samples = sorted(list(dfs_by_sample.keys()))
        series_by_sample = {sample: dfs_by_sample[sample][value_column].rename(sample) for sample in samples}
        if renamed_key_column is None:
            renamed_key_column = key_column
        combined_df = pd.DataFrame(series_by_sample).rename_axis(renamed_key_column).reset_index()
        combined_df.to_csv(output_filepath, sep = sep, index = False)
        return combined_df
    else:
        shell(f'touch {output_filepath}')

# Include different quantification methods
include: 'gene_quantifications_kallisto.smk'
include: 'gene_quantifications_homer.smk'

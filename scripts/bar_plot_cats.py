import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


pd.options.mode.chained_assignment = None

sns.set_theme(style='whitegrid', rc={
    'figure.figsize': (6.4, 6.4),  # figure size in inches
    'axes.labelsize': 14,           # font size of the x and y labels
    'xtick.labelsize': 12,          # font size of the tick labels
    'ytick.labelsize': 12,          # font size of the tick labels
    'legend.fontsize': 12,          # font size of the legend
})

  # set the resolution to 300 DPI

def parse_line(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
    read_length = int(fields[1])
    ref_length = int(fields[6])
    mapping_location = fields[5]
    mapq = fields[11]
    AS_score = None
    NM_score = None
    for field in fields[12:]:
        if field.startswith('NM'):
            NM_score = field.split(':')[2]
        if field.startswith('AS'):
            AS_score = field.split(':')[2]
            break
    return read_id, read_length, ref_length, mapping_location, mapq, AS_score, NM_score

def paf_to_dataframe(file):
    data = [parse_line(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id','read_length','ref_length', 'mapping_location', 'mapq', 'AS_score', 'NM_score'])
    # return a test df
    #
    # df = df.head(1000)
    return df

def filter_max_score_keep_dups(df):
    # Remove all rows with * mapping location
    df = df[df['mapping_location'] != '*']
    # Get rows where 'AS_score' is maximum for each 'read_id' 
    idxmax = df.groupby('read_id')['AS_score'].transform('max') == df['AS_score']
    df_max = df.loc[idxmax]

    # where 'read_id', 'mapping location' and 'AS_score' are duplicated select one of them
    non_unique_max = df_max.duplicated(subset=['read_id', 'mapping_location', 'AS_score'])

    # Return the filtered DataFrame and the number of reads with multiple max scores
    return df_max

def get_multimapping_locations(df):
    df = df[df['mapping_location'] != '*']
    df['mapq'] = df['mapq'].astype(int)
    df = df[df['mapq'] == 0]
    df['AS_score'] = df['AS_score'].astype(float)
    # Get the 'read_id's to remove
    df_max_filter, num_multi_max = filter_max_score(df)
    # Get the mapping locations of the multi-mapped reads
    multi_mapping_locations = df_max_filter['mapping_location'].unique()
    return multi_mapping_locations

def add_haplotype_info(df, output_prefix):
    # Add a column with the haplotype information
    df['mapping_location'] = df['mapping_location'].astype(str).str.strip()
    if 'RIL' in output_prefix:
        df['haplotype'] = df['mapping_location'].str.split('_').str[0]
    if 'Atlantic' in output_prefix:
        pattern = r'(\dG)'
        df['haplotype'] = df['mapping_location'].str.extract(pattern)
    if 'Orang' in output_prefix:
        df['haplotype'] = df['mapping_location'].str.split('_').str[0]
    return df

def add_gene_info(df, output_prefix, gene_list):
    if 'RIL' in output_prefix:
        df['gene'] = df['mapping_location'].str.split('_').apply(lambda x: '_'.join(x[1:]))
    if 'Atlantic' in output_prefix:
        pattern = r'(Synt_\d+)_.*x4'
        df['gene'] = df['mapping_location'].str.extract(pattern)
    if 'Orang' in output_prefix:
        df['gene'] = df['mapping_location'].str.split('_').apply(lambda x: '_'.join(x[1:]))
        df = filter_genes(df, gene_list)
    return df

def pivot_table(df, output_prefix, gene_list):
    df = add_haplotype_info(df, output_prefix)
    df = add_gene_info(df, output_prefix, gene_list)
    # sort on mapping location
    df = df.sort_values(by=['read_id', 'mapping_location'])
    return df.pivot_table(index=['read_id'], values=['mapping_location','mapq', 'AS_score', 'haplotype', 'gene', 'read_length', 'ref_length'], aggfunc=lambda x: ','.join(str(v) for v in x)) # ['sum', 'count']

def check_if_comma(value):
    return ',' not in value

def check_same_values(value):
    # Return False if the value is not contained
    if ',' not in value:
        return False
    else:
        parts = value.split(',')
        # Return False if both values are 'nan'
        if parts[0].strip().lower() == 'nan' and parts[1].strip().lower() == 'nan':
            return False
        return len(set(parts)) == 1

     
def remove_duplicates(s):
    # Convert the string to a set of unique numbers
    numbers = set(float(num_str) for num_str in s.split(','))
    
    # Calculate the mean of the numbers
    mean = np.mean(list(numbers))
    
    # Return the mean as a string
    return int(mean)

def add_mapping_category(pivoted_df):
    pivoted_df['mapping_category'] = 'unclassified'
 
    pivoted_df.loc[pivoted_df['mapping_location'].apply(check_if_comma), 'mapping_category'] = pivoted_df['haplotype'].apply(lambda x: f'ASE_{x}')

    # if mapping location contains a comma, it is a multimapping
    # get the mapping location prefixes
    pivoted_df.loc[~pivoted_df['mapping_location'].apply(check_if_comma), 'mapping_category'] = 'multimapping'

    pivoted_df.loc[pivoted_df['haplotype'].apply(check_same_values), 'mapping_category'] = pivoted_df['haplotype'].str.split(',').str[0].apply(lambda x: f'multimapping_within_{x}')

    # if mapping location contains a comma, it is a multimapping
    # get the mapping location prefixes
    pivoted_df.loc[pivoted_df['gene'].apply(check_same_values), 'mapping_category'] = 'multimapping_same_gene'

    print(pivoted_df)

    return pivoted_df


def add_mapq_flag(pivoted_df):
    pivoted_df['mapq_flag'] = 'NaN'
 
    # pivoted_df.loc[~pivoted_df['mapping_location_seperate'].apply(check_if_comma) & pivoted_df['mapq_seperate'].str.contains("0"), 'mapq_flag'] = 'mapq_0_seperate'
    pivoted_df['mapq'] = pivoted_df['mapq'].apply(remove_duplicates)
    pivoted_df['mapq'] = pivoted_df['mapq'].astype(int)
    pivoted_df.loc[pivoted_df['mapping_location'].apply(check_if_comma) & (pivoted_df['mapq'] < 6), 'mapq_flag'] = 'mapq_lt6'

    # pivoted_df.loc[condition_seperate | condition_competetive, 'mapq_flag'] = 'mapq_0_both'
    return pivoted_df

def make_category_barplot(df, output_prefix, gene_list):
    df = df['mapping_category'].value_counts().reset_index()
    plt.figure(figsize=(10, 6))
    sns.barplot(x='mapping_category', y='count', data=df, ci=None)
    plt.title('Counts in each category')
    plt.xlabel('Mapping Category')
    plt.ylabel('Count')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    plt.tight_layout()

    # Save the plot
    plt.savefig(f'{output_prefix}_categories_barplot.png', dpi=300, bbox_inches='tight')

def get_mapping_cat_per_gene(df):
    # Split the 'mapping_location' column on commas to create a list
    df['gene'] = df['gene'].str.split(',')

    # Use explode to create a new row for each mapping location
    unpivoted_df = df.explode('gene')

    # Group by 'mapping_location' and 'mapping_category' and count the number of reads
    grouped_df = unpivoted_df.groupby(['gene', 'mapping_category']).size().reset_index(name='count')

    #print(grouped_df[grouped_df['mapping_location'].isin(['12272_FBgn0000042', 'w1118_FBgn0000042'])])
    print(grouped_df)
    return grouped_df

def filter_genes(df, gene_list):
    return df[df['gene'].isin(gene_list)]

def calculate_mean(s):
    return np.mean([int(i) for i in s.split(',')])



def write_multi_mapping_file(df, output_prefix):
    # Get the mapping locations of the multi-mapped reads
    #df_multi = df[~df['mapping_category'].str.contains('ASE')]
    df_multi = df
    # group by 'mapping_location' and count the number of reads
    df_multi['read_length'] = df_multi['read_length'].apply(calculate_mean)
    #df_multi['ref_length'] = df_multi['ref_length'].apply(calculate_mean)
    df_multi = df_multi.groupby('mapping_location').agg({'read_length': 'mean', 'ref_length': 'first', 'mapping_category': 'size'}).reset_index()
    df_multi.rename(columns={'mapping_category': 'count', 'read_length': 'read_length_mean'}, inplace=True)

    print(df_multi)
    # save to file
    df_multi.to_csv(f'{output_prefix}_gene_counts.tsv', sep='\t', index=False)
    

def get_haplotype_props(df, output_prefix, gene_list):
    # pivot table to get the counts of each haplotype
    # Calculate the total count for each mapping location
    # Filter df to onlty include unique mapping locations
    df = df[df['mapping_category'].str.contains('ASE')]
    # remove rows with nan values in gene column
    df = df.dropna(subset=['gene'])
    if 'Orang' in output_prefix:
        df = filter_genes(df, gene_list)
    if 'Atlantic' in output_prefix:
        df = df[(df['gene'].str.contains('Synt'))]
    total_counts = df.groupby('gene')['count'].sum()

    # Calculate the proportion of each category for each mapping location
    df['proportion'] = df.apply(lambda row: row['count'] / total_counts[row['gene']], axis=1)
    df = df[df['gene'].map(total_counts) >= 10]
   
    print(df)
    # Sort the DataFrame by 'proportion' in descending order
    df = df.sort_values(by='proportion', ascending=False)
    # Save the DataFrame to a file
    df.to_csv(f'{output_prefix}_haplotype_props.tsv', sep='\t', index=False)
    return df

def filter_mapq(df):
    # Filter the DataFrame on mapq == 0
    return df[df['mapq_flag'] != "mapq_lt6"]

def barplot(df, output_prefix):
    plt.figure(figsize=(10, 6))
    sns.barplot(x='mapping_category', y='count', data=df, ci=None)
    plt.title('Counts in each category')
    plt.xlabel('Mapping Category')
    plt.ylabel('Count')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
    plt.tight_layout()

    # Save the plot
    plt.savefig(f'{output_prefix}_barplot.png', dpi=300, bbox_inches='tight')

def histplot(df, output_prefix):
    # PIVOT DF
    print(df)
    # Drop row with nan in gene colum
    # Set the index back to 'gene'
    # Drop row with nan in gene column
    df = df.dropna(subset=['gene'])
    print(df)
    print(df.columns)
    # Set the index back to 'gene'
    #
    # df = df.set_index('gene')
    df  = df.pivot(index='gene', columns='mapping_category', values='proportion')


    print(df)
    
    if 'RIL' in output_prefix:
        colum_to_plot = 'ASE_12272'
    if 'Atlantic' in output_prefix:
        colum_to_plot = 'ASE_1G'
    if 'Orang' in output_prefix:
        colum_to_plot = 'ASE_hap1'

    plt.figure(figsize=(10, 6))
    
    sns.histplot(data=df, x= colum_to_plot, bins = 100, color = 'blue')
    #plt.title('Histogram of read counts per category')
    plt.xlabel(f'Proportions of {colum_to_plot}')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_histplot.pdf', dpi=300)

def load_data(paf):
    if paf.endswith('.paf'):
        return paf_to_dataframe(paf)
    else:
        return pd.read_csv(paf, sep='\t')

def filter_data(df):
    df = df[df['read_length'] >= 200]
    # get test data
    #df = df.head(10000)
    return filter_max_score_keep_dups(df)

def process_data(df, output_prefix, gene_list):
    pivot = pivot_table(df, output_prefix, gene_list)
    print(pivot)
    pivot_with_cat = add_mapping_category(pivot)
    make_category_barplot(pivot_with_cat, output_prefix, gene_list)
    gene_counts = get_mapping_cat_per_gene(pivot_with_cat)
    write_multi_mapping_file(pivot_with_cat, output_prefix)
    return gene_counts

def main(paf, output_prefix):
    df = load_data(paf)
    gene_list = list(pd.read_csv("/blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/hap_genomes/AG06213_PAB.hap1.hap2.shared_genes.tsv", index_col=0).index)
    filtered_df = filter_data(df)
    
    gene_counts = process_data(filtered_df, output_prefix, gene_list)
    print(gene_counts)
    gene_counts.fillna(0, inplace=True)
    gene_props = get_haplotype_props(gene_counts, output_prefix, gene_list)

    

    print(gene_props.loc[gene_props['proportion'].nlargest(10).index])

    histplot(gene_props, output_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process two PAF files.')
    parser.add_argument('--paf', type=str, help='Path to the first seperate PAF file')
    parser.add_argument('--output_prefix', type=str, help='Prefix output file')
    args = parser.parse_args()
    main(args.paf,  args.output_prefix)


#python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/bar_plot_cats.py --paf test_data/align_seperate_dm12272_01h_rep1.paf  --output_prefix dm12272_01h_rep1_test_RIL





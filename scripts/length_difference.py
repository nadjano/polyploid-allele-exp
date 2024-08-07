import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import gffutils

pd.options.mode.chained_assignment = None

sns.set_theme(style='ticks', rc={
    'figure.figsize': (8, 4),  # figure size in inches
    'axes.labelsize': 16,           # font size of the x and y labels
    'xtick.labelsize': 16,          # font size of the tick labels
    'ytick.labelsize': 16,          # font size of the tick labels
    'legend.fontsize': 16,          # font size of the legend
    "axes.spines.right": False, 
    "axes.spines.top": False})


  # set the resolution to 300 DPI

  # set the resolution to 300 DPI

def parse_line(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
    read_length = int(fields[1])
    ref_length = int(fields[10])
    mapping_location = fields[5]
    mapq = fields[11]
    AS_score = None
    NM_score = None
    for field in fields[12:]:
        if field.startswith('NM'):
            NM_score = field.split(':')[2]
        if field.startswith('ms'):
            AS_score = field.split(':')[2]
            break
    return read_id, read_length, ref_length, mapping_location, mapq, AS_score, NM_score

def paf_to_dataframe(file):
    data = [parse_line(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id','read_length','ref_length', 'mapping_location', 'mapq', 'AS_score', 'NM_score'])
    # return a test df
    #df = df.head(1000)
    return df

def get_transcript_lengths(gff_file, output_prefix):
    # Create a database from the GFF file
    db = gffutils.create_db(gff_file, dbfn='gff.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)

    # Initialize an empty dictionary to store transcript lengths
    transcript_lengths = {}

    # Iterate over all features of type 'exon' in the GFF file
    for exon in db.features_of_type('exon'):

        parent_id = exon.attributes['Parent'][0]
        print(parent_id )
        length = exon.end - exon.start + 1
        print(length)
        if parent_id not in transcript_lengths:
            transcript_lengths[parent_id] = 0

        transcript_lengths[parent_id] += length
    print(transcript_lengths)

    # Filter based on output_prefix if necessary
    if 'Atlantic' in output_prefix:
        transcript_lengths = {k: v for k, v in transcript_lengths.items() if 'Synt' in k and 'x4' in k}
    elif 'Orang' in output_prefix:
        pass  # No additional filtering needed for 'Orang'

    # Create a DataFrame from the dictionary
    df = pd.DataFrame.from_dict(transcript_lengths, orient='index', columns=['ref_length'])
    return df

def pivot_length_table(df, output_prefix):
    df['mapping_location'] = df.index
    df = add_gene_info(df, output_prefix)
    df = add_haplotype_info(df, output_prefix)
    
    df = df.pivot(index='gene', columns='haplotype', values=['ref_length'])
    df.columns = ['_'.join(col) for col in df.columns.values]
    # Make values in the 'ref_length' columns integers
    print(df)
    # Drop rows with NaN values
    df.dropna(inplace=True)
    df.fillna(0, inplace=True)
    df = df.astype(int)
    print(df)
    return df

def get_length_proportions(df, output_prefix):
    # Add colum with synt id and haplotype
    df['mapping_location'] = df.index
    df = add_gene_info(df, output_prefix)
    df = add_haplotype_info(df, output_prefix)
    total_length = df.groupby('gene')['ref_length'].sum()

    df['length_proportion'] = df.apply(lambda row: row['ref_length'] / total_length[row['gene']], axis=1)

    df = df.pivot(index='gene', columns='mapping_category', values=['length_proportion'])
    df.columns = ['_'.join(col) for col in df.columns.values]
    
    print(df)
    return df

def check_length_values(row, percent):
    min_value = row.min()
    max_value = row.max()

    # Calculate the 1% range
    ten_percent = min_value * percent * 0.01

    # Check if all values are within 1% of each other
    return max_value >= min_value + ten_percent

def get_haplotype_with_longest_annotation(row):
    if row['ref_length_hap1'] == row['ref_length_hap2']:
        return 'equal_lengths'
    elif row['ref_length_hap1'] > row['ref_length_hap2']:
        return 'ref_length_hap1'
    else:
        return 'ref_length_hap2'

def add_longest_transcript(df, output_prefix):
    if 'Atlantic' in output_prefix:
        df['haplotype_with_longest_annotation'] = df[['ref_length_haplotype1', 'ref_length_haplotype2', 'ref_length_haplotype3', 'ref_length_haplotype4']].idxmax(axis=1)
        mask = (df[['ref_length_haplotype1', 'ref_length_haplotype2', 'ref_length_haplotype3', 'ref_length_haplotype4']].nunique(axis=1) == 1)
        df.loc[mask, 'haplotype_with_longest_annotation'] = 'equal_lengths'
    if 'Orang' in output_prefix:
        df['haplotype_with_longest_annotation'] = df.apply(get_haplotype_with_longest_annotation, axis=1)
    # Remove the 'ref_length_unique_' prefix
    df['haplotype_with_longest_annotation'] = df['haplotype_with_longest_annotation'].str.replace('ref_length_', '')
    #df.loc[df['length_category'] == 'less_1%_difference', 'haplotype_with_longest_annotation'] = 'lengths_within_5%'
    return df

def add_length_category(df, output_prefix):
    # Select only the relevant columns
    print(df)
    if 'Atlantic' in output_prefix:
        df_subset = df[['ref_length_haplotype1', 'ref_length_haplotype2', 'ref_length_haplotype3', 'ref_length_haplotype4']]
    if 'Orang' in output_prefix:
        df_subset = df[['ref_length_hap1', 'ref_length_hap2']]
    df_subset = df_subset.astype(int)
    print(df_subset)

    # print max and min values of each row
    print(df_subset.apply(lambda x: x.max() - x.min(), axis=1))
    df['length_category'] = 'unclassified'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 0), axis=1), 'length_category'] = 'less_5%_difference'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 5), axis=1), 'length_category'] = 'more_5%_difference'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 10), axis=1), 'length_category'] = 'more_10%_difference'
    df.loc[df_subset.apply(lambda x: check_length_values(x, 20), axis=1), 'length_category'] = 'more_20%_difference'

    print(df['length_category'])
    return df

def add_synt_id(df):
    df['gene'] = df.index.to_series().apply(lambda x: x.split('_chr')[0])
    pattern = r'(\dG)'
    df['haplotype'] = 'ASE_' + df.index.to_series().str.extract(pattern)[0]
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

def add_haplotype_info(df, output_prefix):
    # Add a column with the haplotype information
    df['mapping_location'] = df['mapping_location'].astype(str).str.strip()
    if 'RIL' in output_prefix:
        df['haplotype'] = df['mapping_location'].str.split('_').str[0]
    if 'Atlantic' in output_prefix:
        pattern = r'(\dG)'
        df['haplotype'] = df['mapping_location'].str.extract(pattern)
        df['haplotype'] = df['haplotype'].replace({'0G': 'unphased', '1G': 'haplotype1', '2G': 'haplotype2', '3G': 'haplotype3', '4G': 'haplotype4'})
        print(df['haplotype'])
    if 'Orang' in output_prefix:
        df['haplotype'] = df['mapping_location'].str.split('_').str[0]
    return df

def add_gene_info(df, output_prefix):
    if 'RIL' in output_prefix:
        df['gene'] = df['mapping_location'].str.split('_').apply(lambda x: '_'.join(x[1:]))
    if 'Atlantic' in output_prefix:
        pattern = r'(Synt_\d+)'
        df['gene'] = (df['mapping_location'].str.extract(pattern))
        print(df['gene'].unique())
    if 'Orang' in output_prefix:
        df['gene'] = df['mapping_location'].str.split('_').apply(lambda x: '_'.join(x[1:]))
        print(df['gene'])
    return df

def pivot_table(df, output_prefix):
    df = add_gene_info(df, output_prefix)
    df = add_haplotype_info(df, output_prefix)
    print(df)
    return df.pivot_table(index=['read_id'], values=['mapping_location','mapq', 'AS_score', 'haplotype', 'gene', 'read_length'], aggfunc=lambda x: ','.join(str(v) for v in x)) # ['sum', 'count']

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
    pivoted_df.loc[pivoted_df['gene'].apply(check_same_values), 'mapping_category'] = 'multimapping_same_gene'

    print(pivoted_df[['mapping_location', 'mapping_category']])

    return pivoted_df

def get_mapping_cat_per_gene(df):
    # Split the 'mapping_location' column on commas to create a list
    df['gene'] = df['gene'].str.split(',')

    # Use explode to create a new row for each mapping location
    unpivoted_df = df.explode('gene')

    print('unpivoted_df')
    print(unpivoted_df)
    # Group by 'mapping_location' and 'mapping_category' and count the number of reads
    grouped_df = unpivoted_df.groupby(['gene', 'mapping_category']).size().reset_index(name='count')

    return grouped_df

def filter_genes(df, gene_list):
    return df[df['gene'].isin(gene_list)]

def get_haplotype_props(df, output_prefix, gene_list):
    # pivot table to get the counts of each haplotype
    # Calculate the total count for each mapping location
    # Filter df to onlty include unique mapping locations
    
    if 'Orang' in output_prefix:
        df = filter_genes(df, gene_list)
    total_counts = df.groupby('gene')['count'].sum()

    # Calculate the proportion of each category for each mapping location
    df['ratio'] = df.apply(lambda row: row['count'] / total_counts[row['gene']], axis=1)
    
    df = df[df['gene'].map(total_counts) >= 10]

    return df

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

def pivot_genes(df):
    df  = df.pivot(index='gene', columns='mapping_category', values=['ratio'])
    df.columns = ['_'.join(col) for col in df.columns.values]
    # convert coutn_proportion 
    print(df)
    return df

def box_plot(df, output_prefix):
    # Create a scatter plot of 'prop_ASE_read_G1' vs 'prop_ASE_read_G2'
    #plt.figure(figsize=(10, 10))
    ax = sns.boxplot(data=df, x="length_category", y='ratio', hue='mapping_category')

    # Add number of observations
    for i in range(len(df['length_category'].unique())):
        nobs = df['length_category'].value_counts().values[i]
        ax.text(i, 1.2, f"N: {nobs}", horizontalalignment='center', size='large', color='black', weight='semibold')

    # turn x labels
    plt.xticks(rotation=90)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_boxplot.pdf', dpi=300)
    plt.close()

def boxplot1(df, output_prefix):
    plt.figure(figsize=(4, 4))
    # sort by haplotype with longest annotation
    df = df.sort_values(by='haplotype_with_longest_annotation')
    # Make a boxplot
    ax1 = sns.boxplot(data=df, x="haplotype_with_longest_annotation", y='ratio', hue='mapping_category')

    # Add number of observations
    # for i in range(len(df['haplotype_with_longest_annotation'].unique())):
    #     nobs = df['haplotype_with_longest_annotation'].value_counts().values[i]
    #     ax1.text(i, 1.2, f"N: {nobs}", horizontalalignment='center', size='large', color='black', weight='semibold')
    # turn x labels
    plt.xticks(rotation=90)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()


    plt.savefig(f'{output_prefix}_boxplot2.pdf', dpi=300, bbox_inches='tight', transparent=True)
    plt.close()

def get_gff_file(output_prefix):
    if 'Atlantic' in output_prefix:
        return '/blue/mcintyre/share/potato_ASE/spuddb_reference_data/ATL_v3.hc_gene_models_syntIDs.repr.gff3'
    if 'Orang' in output_prefix:
        return '/blue/mcintyre/share/potato_ASE/orang_utan_references/AG06213_PAB.hap1.hap2.refseq.liftoff.v2_only_gene.gff'

def load_data(paf, output_prefix):
    if paf.endswith('.paf'):
        df = paf_to_dataframe(paf) 
    else:
        df = pd.read_csv(paf, sep='\t')
    return df

def filter_data(df, output_prefix):
    if 'Atlantic' in output_prefix:
        df = df[df['mapping_location'].str.contains('Synt') & df['mapping_location'].str.contains('x4')]
    filtered_df = filter_max_score_keep_dups(df)
    return filtered_df

def process_data(df, output_prefix):
    pivot = pivot_table(df, output_prefix)
    pivot_with_cat = add_mapping_category(pivot)
    gene_counts = get_mapping_cat_per_gene(pivot_with_cat)
    gene_counts_ASE = gene_counts[gene_counts['mapping_category'].str.contains('ASE')]
    return gene_counts_ASE

def main(paf, output_prefix):
    gff_file = get_gff_file(output_prefix)
    transcript_lengths = get_transcript_lengths(gff_file, output_prefix)
    transcript_lengths = pivot_length_table(transcript_lengths, output_prefix)
    transcript_lengths_with_cat = add_length_category(transcript_lengths, output_prefix)
    transcript_lengths_with_cat = add_longest_transcript(transcript_lengths_with_cat, output_prefix)

    df = load_data(paf, output_prefix)
    gene_list = list(pd.read_csv("/blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/hap_genomes/AG06213_PAB.hap1.hap2.shared_genes.tsv", index_col=0).index)
    filtered_df = filter_data(df, output_prefix)
    gene_counts_ASE = process_data(filtered_df, output_prefix)
    ASE_gene_props = get_haplotype_props(gene_counts_ASE, output_prefix, gene_list)
    
    pivoted_genes = pivot_genes(ASE_gene_props)
    count_columns = [col for col in pivoted_genes.columns if 'ratio' in col]
    pivoted_genes = pivoted_genes.dropna(subset=count_columns, thresh=len(count_columns) - 3)
    pivoted_genes = pivoted_genes.reset_index()
    pivoted_genes = pivoted_genes.melt(id_vars=['gene'], value_vars=count_columns, var_name='mapping_category', value_name='ratio')

    joined = pd.merge(pivoted_genes, transcript_lengths, on='gene', how='left')
    print(joined)
    box_plot(joined, output_prefix)
    boxplot1(joined, output_prefix)
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process two PAF files.')
    parser.add_argument('--paf', type=str, help='Path to the PAF file')
    parser.add_argument('--output_prefix', type=str, help='Prefix output file')
    args = parser.parse_args()
    main(args.paf,  args.output_prefix)


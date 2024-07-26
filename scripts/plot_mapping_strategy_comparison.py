import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


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
        if field.startswith('AS'):
            AS_score = field.split(':')[2]
            break
    return read_id, read_length, ref_length, mapping_location, mapq, AS_score, NM_score

def paf_to_dataframe(file):
    data = [parse_line(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id','read_length','ref_length', 'mapping_location', 'mapq', 'AS_score', 'NM_score'])
    # return a test df
    #df = df.head(1000)
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
        # rename the haplotypes 0G to unphased, 1G to hap1, 2G to hap2, 3G to hap3, 4G to hap4
        df['haplotype'] = df['haplotype'].replace({'0G': 'unphased', '1G': 'hap1', '2G': 'hap2', '3G': 'hap3', '4G': 'hap4'})
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

def add_longest_transcript(df, output_prefix):
    if 'Atlantic' in output_prefix:
        df['haplotype_with_longest_annotation'] = df[['ref_length_hap1', 'ref_length_hap2', 'ref_length_hap3', 'ref_length_hap4']].idxmax(axis=1)
    if 'Orang' in output_prefix:
        df['haplotype_with_longest_annotation'] = df[['ref_length_hap1', 'ref_length_hap2']].idxmax(axis=1)
    # Remove the 'ref_length_unique_' prefix
    df['haplotype_with_longest_annotation'] = df['haplotype_with_longest_annotation'].str.replace('ref_length_', '')
    df.loc[df['length_category'] == 'less_5%_difference', 'haplotype_with_longest_annotation'] = 'lengths_within_5%'
    return df

def pivot_table(df, output_prefix):
    df = add_haplotype_info(df, output_prefix)
    df = add_gene_info(df, output_prefix)
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

def join_dfs(dataframes, type, output_prefix):
        # Ensure 'mapping_location' column in both dataframes is not a list
    
    merged = pd.merge(dataframes[0], dataframes[1], on=['mapping_location', 'mapping_category'], suffixes=('_seperate', '_competetive'), how= type)

    # print the rows where the count is different
    diff = merged[merged['count_seperate'] != merged['count_competetive']]
    # Save to file
    diff.to_csv(f'{output_prefix}_diff.tsv', sep='\t', index=False)

    return merged
     
def remove_duplicates(s):
    # Convert the string to a set of unique numbers
    numbers = set(float(num_str) for num_str in s.split(','))
    
    # Calculate the mean of the numbers
    mean = np.mean(list(numbers))
    
    # Return the mean as a string
    return int(mean)

def add_mapping_category(pivoted_df):
    pivoted_df['mapping_category'] = 'multimapping_multi_genes'
 
    pivoted_df.loc[pivoted_df['mapping_location'].apply(check_if_comma), 'mapping_category'] = pivoted_df['haplotype'].apply(lambda x: f'unique_{x}')

    # if mapping location contains a comma, it is a multimapping
    # get the mapping location prefixes
    pivoted_df.loc[pivoted_df['gene'].apply(check_same_values), 'mapping_category'] = 'multimapping_same_gene'

    print(pivoted_df[['mapping_location', 'mapping_category']])

    #cat_counts = get_mapping_cat_per_gene(pivoted_df)

    pivoted_df = add_mapq_flag(pivoted_df)
    print(pivoted_df[['mapping_location', 'mapq', 'mapq_flag']])

    return pivoted_df


def add_mapq_flag(pivoted_df):
    pivoted_df['mapq_flag'] = 'NaN'
 
    # pivoted_df.loc[~pivoted_df['mapping_location_seperate'].apply(check_if_comma) & pivoted_df['mapq_seperate'].str.contains("0"), 'mapq_flag'] = 'mapq_0_seperate'
    pivoted_df['mapq'] = pivoted_df['mapq'].apply(remove_duplicates)
    pivoted_df['mapq'] = pivoted_df['mapq'].astype(int)
    pivoted_df.loc[pivoted_df['mapping_location'].apply(check_if_comma) & (pivoted_df['mapq'] < 6), 'mapq_flag'] = 'mapq_lt6'

    # pivoted_df.loc[condition_seperate | condition_competetive, 'mapq_flag'] = 'mapq_0_both'
    return pivoted_df


def get_mapping_cat_per_gene(df):
    # Split the 'mapping_location' column on commas to create a list
    df['mapping_location'] = df['mapping_location'].str.split(',')

    # Use explode to create a new row for each mapping location
    unpivoted_df = df.explode('mapping_location')

    # Group by 'mapping_location' and 'mapping_category' and count the number of reads
    grouped_df = unpivoted_df.groupby(['mapping_location', 'mapping_category']).size().reset_index(name='count')

    #print(grouped_df[grouped_df['mapping_location'].isin(['12272_FBgn0000042', 'w1118_FBgn0000042'])])

    return grouped_df


def get_read_with_different_mappings(dataframes, output_prefix):
    merged = pd.merge(dataframes[0], dataframes[1], on='read_id', suffixes=('_seperate', '_competetive'))
 
    # print the rows where the count is different
    diffrent_count = merged[merged['mapping_location_seperate'] != merged['mapping_location_competetive']]
    # Save to file
    diffrent_count.to_csv(f'{output_prefix}_read_diff.tsv', sep='\t', index=True)

    diff_group = diffrent_count.groupby(['mapping_category_seperate','mapping_category_competetive']).size().reset_index(name='counts')
    #print(diff_group)

    # save to file
    diff_group.to_csv(f'{output_prefix}_grouped_categories.tsv', sep='\t', index=False)
    return merged

def create_scatter_plot(df, output_prefix):
    # Create a color dictionary to map each unique mapping category to a color
    # get the unique mapping categories
    mapping_categories = df['mapping_category'].unique()
    df['mapping_category'] = df['mapping_category'].replace({'multimapping_multi_genes': 'Multimapping, multi genes', 'multimapping_same_gene': 'Multimapping, same gene', 'unique_12272': 'Unique mapping, haplotype 1', 'unique_w1118': 'Unique mapping, haplotype 2', 'unique_hap1': 'Unique mapping, haplotype 1', 'unique_hap2': 'Unique mapping, haplotype 2', 'unique_hap3': 'Unique mapping, haplotype 3', 'unique_hap4': 'Unique mapping, haplotype 4', 'unique_unphased': 'Unique mapping, unphased'})
    #print(mapping_categories)
    # Create a dictionary to map each unique mapping category to a color and a shape
    if "RIL" in output_prefix:
        my_markers = {'Multimapping, multi genes': "X", "Multimapping, same gene": "s", 'Unique mapping, haplotype 1': "o", 'Unique mapping, haplotype 2': "o"}
        
        my_palette = {'Multimapping, multi genes': "grey", 'Multimapping, same gene': "black", 'Unique mapping, haplotype 1': "red", 'Unique mapping, haplotype 2': "blue"}

    if "Orang" in output_prefix:
        my_markers = {'Multimapping, multi genes': "X", 'Multimapping, same gene': "s", 'Unique mapping, haplotype 1': "o", 'Unique mapping, haplotype 2': "o"}
        
        my_palette = {'Multimapping, multi genes': "grey", 'Multimapping, same gene': "black", 'Unique mapping, haplotype 1': "red", 'Unique mapping, haplotype 2': "blue"}

    if "Atlantic" in output_prefix:
        my_markers = {'Multimapping, multi genes': "X", 'Multimapping, same gene': "s", 'Unique mapping, haplotype 1': "o", 'Unique mapping, haplotype 2': "o", 'Unique mapping, haplotype 3': "o", 'Unique mapping, haplotype 4': "o",  'Unique mapping, unphased': "o"}
        
        my_palette = {'Multimapping, multi genes': "grey",'Multimapping, same gene': "black", 'Unique mapping, haplotype 1': "red", 'Unique mapping, haplotype 2': "blue", 'Unique mapping, haplotype 3': "green", 'Unique mapping, haplotype 4': "orange", 'Unique mapping, unphased': "purple"}

    # Get unique categories
    categories = df['mapping_category'].unique()

    # Create a PdfPages object
    # with PdfPages(f'{output_prefix}_scatterplot.pdf') as pdf:
    #     # Create a scatterplot for each category
    #     for category in categories:
    #         df_subset = df[df['mapping_category'] == category]
    #         plt.figure(dpi=300)  # Create a new figure
    #         ax = sns.scatterplot(data=df_subset, x="count_seperate", y="count_competetive", hue="mapping_category", style="mapping_category", palette=my_palette, markers=my_markers)
    #         ax.set_xlabel("Seperate mapping count")
    #         ax.set_ylabel("Competetive mapping count")
    #         #ax.set_xlim(0, 30000)
    #         #ax.set_ylim(0, 30000)
    #         plt.title(f'Scatterplot for {category}')  # Set the title of the plot
    #         pdf.savefig()  # Save the current figure to the PDF
    #         plt.clf()  # Clear the current figure
    # Rename the values in mapping category colum
    
    # sort the data on the mapping category
    df = df.sort_values(by='mapping_category')
    scatterplot = sns.scatterplot(data=df, x="count_seperate", y="count_competetive", hue="mapping_category", style="mapping_category", palette= my_palette, markers=my_markers)

    # calculate the correlation between the two columns
    corr = df['count_seperate'].corr(df['count_competetive'])


    # Set axes to log scale
    scatterplot.set_xscale("log")
    scatterplot.set_yscale("log")
    # Rename the legend
    # place legend outside of plot
    scatterplot.legend(title="Mapping category", bbox_to_anchor=(1.05, 1), loc='upper left')

    # make the x and y axis labels the same
    # Add ax labels
    #scatterplot.set_xlabel("Count per gene (Seperate mapping to haplotypes)")
    scatterplot.set_xlabel("")

    #scatterplot.set_ylabel("Count per gene (Competetive mapping to genotype)")
    scatterplot.set_ylabel("")
    # add the correlation to the plot
    scatterplot.text(0.35, 0.25, f'corr: {corr:.2f}', fontsize=16, transform=plt.gcf().transFigure)

    # save the plot
    fig = scatterplot.get_figure()
    fig.tight_layout()
    fig.savefig(f"all_{output_prefix}_scatter_plot.pdf", dpi=300, bbox_inches='tight', transparent=True)

def filter_mapq(df):
    # Filter the DataFrame on mapq == 0
    return df[df['mapq_flag'] != "mapq_lt6"]

def create_plot_and_close(pivot, join_type, output_prefix):
    pivot = [df.copy() for df in pivot]
    gene_counts = [get_mapping_cat_per_gene(df) for df in pivot]
    print(gene_counts)
    merged = join_dfs(gene_counts, join_type, output_prefix)
    print(merged)
    create_scatter_plot(merged, f"{output_prefix}")
    plt.close()


def main(paf1, paf2, output_prefix):
    if paf1.endswith('.paf'):
        dfs = [paf_to_dataframe(paf) for paf in [paf1, paf2]]
    else:
        # read in pafs as tsv files
        dfs = [pd.read_csv(paf, sep='\t') for paf in [paf1, paf2]]
        # set a subset for testing
        #dfs = [df.head(10000) for df in dfs]
 
    # dfs = [paf_to_dataframe(paf) for paf in [paf1, paf2]]
    # Filter dfs for minimum read length
    dfs = [df[df['read_length'] >= 200] for df in dfs]
    # Sort df on mapping_location
    dfs = [df.sort_values(by='mapping_location') for df in dfs]
    pivot_with_cat = []
    gene_counts = []
    
    for df in dfs:
        filtered_df = filter_max_score_keep_dups(df)
        # Sort the df on mapping location
        filtered_df = filtered_df.sort_values(by='mapping_location')
        pivot = pivot_table(filtered_df, output_prefix)
        pivot_with_cat.append(add_mapping_category(pivot))


    # create a scatter plot for unfiltered data
    create_plot_and_close(pivot_with_cat, "outer", output_prefix)

    # Get the reads that have different mappings
    merged_diff = get_read_with_different_mappings(pivot_with_cat, output_prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process two PAF files.')
    parser.add_argument('--paf_seperate', type=str, help='Path to the first seperate PAF file')
    parser.add_argument('--paf_competetive', type=str, help='Path to the first seperate PAF file')
    parser.add_argument('--output_prefix', type=str, help='Prefix output file')
    args = parser.parse_args()
    main(args.paf_seperate, args.paf_competetive, args.output_prefix)


#python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/plot_mapping_strategy_comparison.py --paf_seperate out_ms/paf_tsv_files/RIL_updated_seperate_P.tsv  --paf_competetive out_ms/paf_tsv_files/RIL_updated_competetive_P.tsv --output_prefix RIL_test_size
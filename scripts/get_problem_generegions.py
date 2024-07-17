import argparse
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
import seaborn as sns

# python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/get_problem_generegions.py --pafs align_seperate_SRR14993892_hap1.paf align_seperate_SRR14993892_hap2.paf align_competetive_SRR14993892_allhap.paf 

pd.options.mode.chained_assignment = None

def parse_line_paf(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
    read_length = fields[1]
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
    #return read_id, read_length, mapping_location, mapq, AS_score, NM_score
    return read_id, mapping_location, AS_score

def paf_to_dataframe(file):
    # Only select first 500 lines for testing
    #data = [parse_line_paf(line) for line in open(file, 'r')]
    # Only select first 500 lines for testing
    with open(file, 'r') as f:
        data = [parse_line_paf(line) for line in f if 'x4' in line]
    # Use part of the read name as colum suffix

    df = pd.DataFrame(data, columns=['read_id', 'mapping_location', 'AS_score'])
    read_name = extract_hap_from_filename(file)

    df.columns = [f"{col}_{read_name}" for col in df.columns]
    

    return df

def get_frequencies(df, file):
    read_name = extract_hap_from_filename(file)
    #df = df[1:1000]
    df_counts = df[f"mapping_location_{read_name}"].value_counts().reset_index()
    df_counts.columns = ['mapping_location', f"count_{read_name}"]

    return df_counts

def extract_hap_from_filename(filename):
    # Split the filename into name and extension
    name, ext = os.path.splitext(filename)
    
    # Split the name by underscore and get the last part
    hap = name.split('_')[-1]
    
    return hap

def merge_frequencies(dfs):
    # Join data
    result = dfs[0]

    for df in dfs[1:]:
        result = result.merge(df, on='mapping_location', how='outer')

    # Create a new column 'id' that contains the common part of the IDs
    result['id'] = result['mapping_location'].str.rsplit('chr').str[0]

    # Group by the new column 'id' and sum the values
    result = result.groupby('id').sum().reset_index()
    print(result)
    # Sort the columns
    #result = result.reindex(sorted(result.columns), axis=1)

    # Set 'id' as the index
    result.set_index('id', inplace=True)

    # Save result
    #result.to_csv(f"{out_dir}/{sra_id}_joined_synt.csv")

    return result

def pivot_table(df):
    # Create a new column 'id' that contains the common part of the IDs
    df['id'] = df['mapping_location'].str.rsplit('chr').str[0]
    df['hap'] = df['mapping_location'].str.split('_').str[3]
    # Create a pivot table
    df = df.pivot_table(index='id', columns='hap', values='count_allhap')

    return df

def group_allels_by_mappings(df):
    """Group data based on AS score > 0 pattern."""

    # Group the reads ids based on the pattern of AS scores > 0
    df = df.fillna(0)
    print(df)
    # Create a new column 'id' that contains the common part of the IDs
    df = pivot_table(df)
    print(df)
    cols = [f'{i}G' for i in range(1,5)]
    df['pattern'] = df[cols].apply(lambda row: 'E' if max(row) - min(row) <= 0.2 * max(row) else '|'.join(['1' if x > 1 else '0' for x in row]), axis=1)
    print(df)
    #df.to_csv(f"{out_dir}/{sra_id}_pattern_synt.csv")
    # Group by the 'pattern' column and count the frequencies
    pattern_counts = df.groupby('pattern').size()
    print(pattern_counts)
    return df
    # group_ids = df.gt(0).astype(int).diff().ne(0).cumsum(axis=1).apply(tuple, axis=1)

    # # Get for each group the column names that are True 
    # group_columns = df.gt(0).groupby(group_ids).apply(lambda x: x.columns[x.any()])

    # group_columns_dict = group_columns.apply(list).to_dict()

    # # Now we extract the unique values from the dictionary
    # unique_values = set(tuple(v) for v in group_columns_dict.values())
    # # Write the unique values to a file
    # with open(f"{outdir}/{sra_id}_unique_groups.txt", "w") as f:
    #     for value in unique_values:
    #         f.write(f"{value}\n")


def group_allels_by_patterns(df):

    df = df[df['mapping_location'].str.contains('x4')]
    # Include the genotype and phenotype column
    cols = [f'count_hap{i}' for i in range(1,5)] + [f'count_allhap'] 

    # Create a new column 'pattern' that represents the pattern of 0 and >0 for each row
    df['pattern'] = df[cols].apply(lambda row: 'E' if max(row) - min(row) <= 0.5 * max(row) else '|'.join(['1' if x > 1 else '0' for x in row]), axis=1)

    print(df)
    #df.to_csv(f"{out_dir}/{sra_id}_pattern_synt.csv")
    # Group by the 'pattern' column and count the frequencies
    pattern_counts = df.groupby('pattern').size()
    pattern_counts.columns = ['H1H2H3H4G', 'count']
    pattern_count11111 = df[df['pattern'] == '11111']
    print(pattern_count11111)
    print(pattern_counts)
    # Check for the dataframe how many of the phenotype values are the same
    print(df[df['pattern'] == '0|1|0|0|1'])

    # Save the resultgene
    #pattern_count11111_all_different.to_csv(f"{out_dir}/{sra_id}_pattern_group_countsynt_all_different.csv")

    #pattern_counts.to_csv(f"{out_dir}/{sra_id}_pattern_group_countsynt.csv")


def filter_max_score(df):
    # Get rows where 'AS_score' is maximum for each 'read_id'
    idxmax = df.groupby('read_id_allhap')['AS_score_allhap'].transform('max') == df['AS_score_allhap']
    df_max = df.loc[idxmax]

    # Get rows where 'read_id' and 'AS_score' are duplicated
    non_unique_max = df_max.duplicated(subset=['read_id_allhap', 'AS_score_allhap'], keep=False)
    df_multi = df_max[non_unique_max]

    # Get the unique read ids in the multi-mapped reads
    df_multi_unique = df_multi['read_id_allhap'].unique()

    # Filter out the non-unique max scores
    df_max_filter = df_max[~non_unique_max]

    # Return the filtered DataFrame and the number of reads with multiple max scores
    return df_max_filter

def add_ASE_proportions(df):
    # Replace NaN with 0
    df = df.fillna(0)

    df["proportion_hap1"] = df['1G'] / (df['1G'] + df['2G'] + df['3G'] + df['4G'])
    df["proportion_hap2"] = df['2G'] / (df['1G'] + df['2G'] + df['3G'] + df['4G'])
    df["proportion_hap3"] = df['3G'] / (df['1G'] + df['2G'] + df['3G'] + df['4G'])
    df["proportion_hap4"] = df['4G'] / (df['1G'] + df['2G'] + df['3G'] + df['4G'])

    return df

def adapt_proportions(df):
    # Replace NaN with 0
    df = df.fillna(0)

    # df["proportion_hap1"] =  df["proportion_hap1"] - 0.25
    # df["proportion_hap2"] =  df["proportion_hap2"] - 0.25
    # df["proportion_hap3"] =  df["proportion_hap3"] - 0.25
    # df["proportion_hap4"] =  df["proportion_hap4"] - 0.25

    df["coordinate_x"] =  df["proportion_hap1"] - df["proportion_hap3"]
    df["coordinate_y"] =  df["proportion_hap2"] - df["proportion_hap4"]

    return df

def scatter_plot(df, lable):
    # Plot the proportions of ASE reads per haplotype

    print(df)
    sns.set_theme(style="whitegrid")
  
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x="coordinate_x", y="coordinate_y", data=df, palette='viridis')
    
    # Set labels and title
    plt.xlabel('Haplotypes')
    plt.ylabel('Proportion of ASE reads')

    # Save plot
    plt.savefig(f"scatter_{lable}.png")

def plot_proportions(df):
    # Plot the proportions of ASE reads per haplotype
    df = df.reset_index()
    df = df['id', 'proportion_hap1', 'proportion_hap2', 'proportion_hap3', 'proportion_hap4']
    print(df)
    sns.set_theme(style="whitegrid")
  
    plt.figure(figsize=(10, 6))
    sns.violinplot(x='haplotype', y='value', hue='haplotype', data=df, palette='viridis', inner=None)
    sns.stripplot(x='haplotype', y='value', hue='haplotype', data=df, color='black', size=1)
    
    # Set labels and title
    plt.xlabel('Haplotypes')
    plt.ylabel('Proportion of ASE reads')

    # Save plot
    plt.savefig("proportions.png")

def make_ridge_plot(df):
    # Plot the proportions of ASE reads per haplotype
    df = df.reset_index()
    df = df.melt(id_vars='gene', value_vars=['proportion_hap1', 'proportion_hap2', 'proportion_hap3', 'proportion_hap4'])
    print(df)
  
    plt.figure(figsize=(10, 6))
   
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
    palette = sns.color_palette("Set2", 12)
    g = sns.FacetGrid(df, palette=palette, row="haplotype", hue="haplotype", aspect=9, height=1.2)
    g.map_dataframe(sns.kdeplot, x="value", fill=True, alpha=1)
    g.map_dataframe(sns.kdeplot, x="value", color='black')
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, color='black', fontsize=13,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, "haplotype")
    g.fig.subplots_adjust(hspace=-.5)
    g.set_titles("")
    g.set(yticks=[], xlabel="value")
    g.despine( left=True)
    plt.suptitle('proportions', y=0.98)
    plt.savefig("_ridge_proportions.png")

def main(pafs):
    # Read in the paf files as list
    dfs = [paf_to_dataframe(paf) for paf in pafs]

    df_counts = []

    for df, paf in zip(dfs, pafs):
        df_counts.append(get_frequencies(df, paf))

    # Merge the dataframes
    merged = merge_frequencies(df_counts)
    group_allels_by_patterns(merged)

    df_comp_seondary_groups = group_allels_by_mappings(df_counts[4])

    # Get unique counts for competetive alignent
    df_ASE = filter_max_score(dfs[4])
    # SUM counts per gene
    df_ASE_counts = get_frequencies(df_ASE, pafs[4])

    # Group syntelogs
    df_ASE_counts_grouped = group_allels_by_mappings(df_ASE_counts)

    ASE_with_props = add_ASE_proportions(df_ASE_counts_grouped)

    ASE_with_props = adapt_proportions(ASE_with_props)

    scatter_plot(ASE_with_props,  'unfiltered')

    # Now filter out the syntelogs that have now 1111 pattern 
    df_comp_seondary_groups_1111 = df_comp_seondary_groups[df_comp_seondary_groups['pattern'] == '1|1|1|1']

    # Filter the df_ASE_counts_grouped based on the rows from df_comp_seondary_groups
    ASE_with_props_filtered = ASE_with_props[ASE_with_props.index.isin(df_comp_seondary_groups_1111.index)]

    print(ASE_with_props_filtered)
    scatter_plot(ASE_with_props_filtered, 'filtered')

    df_comp_seondary_groups_E = df_comp_seondary_groups[df_comp_seondary_groups['pattern'] == 'E']

    # Filter the df_ASE_counts_grouped based on the rows from df_comp_seondary_groups
    ASE_with_props_filteredE = ASE_with_props[ASE_with_props.index.isin(df_comp_seondary_groups_E.index)]

    print(ASE_with_props_filtered)
    scatter_plot(ASE_with_props_filteredE, 'filteredE')
    print(ASE_with_props_filteredE)

    # Filter for pattern 1111
    ASE_with_props_filtered12 = ASE_with_props_filtered[ASE_with_props_filtered['pattern'] == '1|1|1|1']
    scatter_plot(ASE_with_props_filtered12, 'filtered12')

    ASE_with_props_filtered2 = ASE_with_props[ASE_with_props['pattern'] == '1|1|1|1']
    scatter_plot(ASE_with_props_filtered2, 'filtered2')



    






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process two PAF files.')
    # Read in the paf files as list
    parser.add_argument('--pafs', type=str, nargs='+', help='Names of the experiments')
    args = parser.parse_args()


    main(args.pafs)
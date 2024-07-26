import pandas as pd
import numpy as np
import os
import sys
import argparse
import re
import matplotlib.pyplot as plt
import seaborn as sns

pd.options.mode.chained_assignment = None

sns.set_theme(style='ticks', rc={
    'figure.figsize': (3, 3),  # figure size in inches
    'axes.labelsize': 14,           # font size of the x and y labels
    'xtick.labelsize': 14,          # font size of the tick labels
    'ytick.labelsize': 14,          # font size of the tick labels
    'legend.fontsize': 14,          # font size of the legend
    "axes.spines.right": False,
    "axes.spines.top": False})

  # set the resolution to 300 DPI


def add_haplotype_info(df, output_prefix):
    # Add a column with the haplotype information
    #df['mapping_location'] = df['mapping_location'].astype(str).str.strip()
    if 'RIL' in output_prefix:
        df['haplotype'] = df['mapping_location'].apply(lambda x: [i.split('_')[0] for i in x])
    if 'Atlantic' in output_prefix:
        pattern = r'(\dG)'
        #df['haplotype'] = df['mapping_location'].str.extract(pattern)
        df['haplotype'] = df['mapping_location'].apply(lambda x: [re.search(pattern, i).group(1) if re.search(pattern, i) else np.nan for i in x])
    if 'Orang' in output_prefix:
        df['haplotype'] = df['mapping_location'].apply(lambda x: [i.split('_')[0] for i in x])
    return df

def add_gene_info(df, output_prefix):
    if 'RIL' in output_prefix:
        df['gene'] = df['mapping_location'].apply(lambda x: [i.split('_',1)[1] for i in x])
    if 'Atlantic' in output_prefix:
        # only get the genes if they are on all 4 haplotypes
        pattern = r'(Synt_\d+)_.*x4'
        #df['gene'] = df['mapping_location'].str.extract(pattern)
        df['gene'] = df['mapping_location'].apply(lambda x: [re.search(pattern, i).group(1) if re.search(pattern, i) else np.nan for i in x])
    if 'Orang' in output_prefix:
        df['gene'] = df['mapping_location'].apply(lambda x: [i.split('_', 1)[1] for i in x])
    return df

def process_gene_counts(gene_counts, output_prefix):
    df = pd.read_csv(gene_counts, sep='\t')
    print(df)
    df['mapping_location'] = df['mapping_location'].str.split(',')
    print(df.explode('mapping_location'))

    df = add_haplotype_info(df, output_prefix)
    df = add_gene_info(df, output_prefix)

    # drop rows with NaN genes
    df = df.dropna(subset=['gene'])
    print(df['gene'])

    # Check if gene appears exactly twice and the genes are the same
    df['gene_count'] = df['gene'].apply(lambda x: len(x) == 2 and x[0] == x[1])

    # If gene appears exactly twice and the genes are the same, set haplotype to 'multi'
    df.loc[df['gene_count'], 'haplotype'] = 'same_gene'

    # Drop the 'gene_count' column as it's no longer needed
    df = df.drop('gene_count', axis=1)

    # Check if there is more than one gene and the genes are not the same
    df['mixed_haplotype'] = df['gene'].apply(lambda x: len(x) > 1 and x[0] != x[1])

    # If there is more than one gene and the genes are not the same, set haplotype to 'mixed'
    df.loc[df['mixed_haplotype'], 'haplotype'] = 'mixed'

    # Drop the 'mixed_haplotype' column as it's no longer needed
    df = df.drop('mixed_haplotype', axis=1)

    print(df.loc[df['haplotype'] == 'mixed'])
    # remove the list in the hapltype colum
    df['haplotype'] = df['haplotype'].apply(lambda x: x[0])
    return df


def get_mappability_per_gene(df, output_prefix):
    unpivoted_df = df.explode('gene')
    print(unpivoted_df)
    # for each gene count the number of rows it has with m haplotype
    unpivoted_df['M_count'] = unpivoted_df['haplotype'].apply(lambda x: x == 'm')
    M_counts = unpivoted_df.groupby('gene')['M_count'].sum()
    print(M_counts)
    pivot = unpivoted_df.pivot_table(index='gene', columns='haplotype', values='count', aggfunc='sum', fill_value=0)
    print(pivot)
    # rename the columns m to multi_genes and s to same_gene
    pivot = pivot.rename(columns={'m': 'multi_genes', 's': 'same_gene'})
    # devide the values in s colum by 2
    pivot['same_gene'] = pivot['same_gene'] / 2
    # print the rows with m > 0
    # sort by m
    pivot = pivot.sort_values(by='multi_genes', ascending=False)
    print(pivot.loc[pivot['multi_genes'] > 5])
    if 'Orang' in output_prefix:
        pivot['mappability'] = (pivot['hap1'] + pivot['hap2'])/ (pivot['multi_genes'] + pivot['same_gene'] +  pivot['hap1'] + pivot['hap2']) 
    if 'RIL' in output_prefix:
        pivot['mappability'] = (pivot['w1118'] + pivot['12272'])/ (pivot['multi_genes'] + pivot['same_gene'] +  pivot['w1118'] + pivot['12272'])
    if 'Atlantic' in output_prefix:
        # rename the columns 1G, 2G, 3G, 4G to hap1, hap2, hap3, hap4
        pivot = pivot.rename(columns={'1G': 'hap1', '2G': 'hap2', '3G': 'hap3', '4G': 'hap4'})
        pivot['mappability'] = (pivot['hap1'] + pivot['hap2'] + pivot['hap3'] + pivot['hap4'])/ (pivot['multi_genes'] + pivot['same_gene'] +  pivot['hap1'] + pivot['hap2'] + pivot['hap3'] + pivot['hap4'])
    print(pivot)
    # LOC100447686
    #save the pivot table to a file
    pivot.to_csv(f'{output_prefix}_mappability_per_gene.csv')
    return pivot

def plot_histogram(df, output_prefix):
    sns.histplot(df['mappability'], bins=50, color='black', )
    #plt.hist(df['mappability'], bins=100)
    # add number of genes to the plot
    plt.text(0.5, 0.5, f'#genes: {len(df)}', fontsize=12, transform=plt.gcf().transFigure)
    plt.xlabel('Mappability')
    plt.ylabel('Number of genes')
    plt.savefig(f'{output_prefix}_mappability_per_gene_min20_counts.pdf', dpi=300, bbox_inches='tight', transparent=True)


def parse_input ():
    parser = argparse.ArgumentParser(description='Get mappability per gene')
    parser.add_argument('-i', '--input', help='Input file with counts per mapping location', required=True)
    parser.add_argument('-o', '--output_prefix', help='Output file with gene names and mappability values', required=True)
    args = parser.parse_args()
    return args

def main():
    args = parse_input()
    gene_counts = process_gene_counts(args.input, args.output_prefix)
    gene_counts = get_mappability_per_gene(gene_counts, args.output_prefix)
    # filter to only include genes with total counts > 20 
    gene_counts = gene_counts.loc[gene_counts.iloc[:, 0:-1].sum(axis=1) > 20]
    plot_histogram(gene_counts, args.output_prefix)

if __name__ == '__main__':
    main()




# python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/get_mappability_per_gene.py -i out_ms/bar_plots/RIL_updated_P_counts_tsv/RIL_updated_P_gene_counts.tsv  -o RIL_ms


#python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/get_mappability_per_gene.py -i out_ms/bar_plots/Orangutan_N200_counts_tsv/Orangutan_N200_gene_counts.tsv  -o Orang-ms

#python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/get_mappability_per_gene.py -i out_ms/bar_plots/Atlantic_N200_counts_tsv/Atlantic_N200_gene_counts.tsv  -o Atlantic-ms
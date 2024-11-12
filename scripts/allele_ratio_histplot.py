import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


pd.options.mode.chained_assignment = None

sns.set_theme(style='ticks', rc={
    'figure.figsize': (4, 4),  # figure size in inches
    'axes.labelsize': 10,           # font size of the x and y labels
    'xtick.labelsize': 10,          # font size of the tick labels
    'ytick.labelsize': 10,          # font size of the tick labels
    'legend.fontsize': 10,          # font size of the legend
    "axes.spines.right": False, 
    "axes.spines.top": False})



def parse_args():
    parser = argparse.ArgumentParser(description='Plot allele ratio histogram')
    parser.add_argument('input', type=str, help='Input file')
    parser.add_argument('output', type=str, help='Output file')
    return parser.parse_args()

# Load input file
def load_input(input_file):
    return pd.read_csv(input_file, sep='\t')

def sum_sample_counts(df):
    # sum the counts for each sample and skip the first column
    df['num_reads'] = df.iloc[:, 1:].sum(axis=1)
    print(df)
    return df


# extract genes of interest
def extract_genes(df):
    df = df[df['tname'].str.contains('4x')]
    print(df)
    return df

def get_multimapping_tags(df):
    # extract the tnames that contain "&"
    df = df[df['tname'].str.contains('&')]
    df['Synt'] = df['tname'].str.extract(r'(Synt_\d+)')
    # sum the number of reads for each Synt
    df = df.groupby('Synt').agg({'num_reads': 'sum'}).reset_index()
    print(df)
    # rename num_reads to multimapping_reads
    df = df.rename(columns={'num_reads': 'multimapping_reads'})
    print(df)
    return df


def make_unique_synt_table(df):
    # extract the tnames that don't contain "&"
    df = df[~df['tname'].str.contains('&')]
    df['Synt'] = df['tname'].str.extract(r'(Synt_\d+)')
    df['haplotype'] = df['tname'].str.extract(r'(hap\d+)')

    # pivot df to get the counts of each haplotype
    df = df.pivot_table(index='Synt', columns='haplotype', values='num_reads', aggfunc='sum').reset_index()
    # drop the rows that have a >0 in the hap0 column
    print(df)
    df = df[df['hap0'].isnull()]
    # drop the hap0 column
    df = df.drop(columns='hap0')
    # remove columns with only 0s except for the Synt column
    # set the Synt column as the index
    df = df.set_index('Synt')
    df = df.loc[:, (df != 0).any(axis=0)]
    return df

def add_mulitmapping_ratio( df, multi):
    # join the two dataframes
    df = df.join(multi.set_index('Synt'), on='Synt')
    # calculate the ratio of multimapping reads vs total reads
    df['multimapping_ratio'] = df['multimapping_reads'] / (df['hap1'] + df['hap2'] + df['hap3'] +  df['hap4']+ df['multimapping_reads'])
    # add a tag if the multimapping ratio is above 0.25
    df['multimapping_tag'] = np.where(df['multimapping_ratio'] > 0.25, 'multimapping', 'unique')
    return df

def get_allele_ratio(df):
    # replace NaN with 0
    df = df.fillna(0)
    df['allele_ratio_hap1'] = df['hap1'] / (df['hap1'] + df['hap2'] + df['hap3'] + df['hap4'])
    return df

def plot_hist(df, output_file):
    # drop allele ratios that are NaN
    df = df.dropna(subset=['allele_ratio_hap1'])
    print(df)
    sns.histplot(data=df, x= f'allele_ratio_hap1', bins = 100, color = 'red', hue='multimapping_tag', multiple='stack')
    plt.xlabel('Allele ratio hap1')
    plt.ylabel('Count')
    plt.xlim(0, 1)
    plt.savefig(output_file)
    plt.close()
    # also save the ratios to a file
    df.to_csv(output_file + '.tsv', sep='\t')


def main():
    args = parse_args()
    df = load_input(args.input)
    df = sum_sample_counts(df)
    multi = get_multimapping_tags(df)
    print(multi)
    df = extract_genes(df)
    df = make_unique_synt_table(df)
    df = get_allele_ratio(df)
    df = add_mulitmapping_ratio(df, multi)

    plot_hist(df, args.output)


if __name__ == '__main__':
    main()
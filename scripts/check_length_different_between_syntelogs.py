import pandas as pd
import gffutils
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import os

# function to read in the gtf file
def get_transcript_lengths(gff_file):
    # Create a database from the GFF file
    if os.path.exists('gff.db'):
        db = gffutils.FeatureDB('gff.db')
    else:
        db = gffutils.create_db(gff_file, dbfn='gff.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True, disable_infer_transcripts=True, disable_infer_genes=True)
    

    # Initialize an empty dictionary to store transcript lengths
    transcript_lengths = {}

    # Iterate over all features of type 'exon' in the GFF file
    for exon in db.features_of_type('exon'):

        parent_id = exon.attributes['transcript_id'][0]
        length = exon.end - exon.start + 1
        if parent_id not in transcript_lengths:
            transcript_lengths[parent_id] = 0

        transcript_lengths[parent_id] += length
    #print(transcript_lengths)



    # Create a DataFrame from the dictionary
    df = pd.DataFrame.from_dict(transcript_lengths, orient='index', columns=['ref_length'])
    print(df)
    return df


def add_haplotype_info(df):
    # Add a column with the haplotype information
    df['mapping_location'] = df['mapping_location'].astype(str).str.strip()
    pattern = r'(hap\d)'
    df['haplotype'] = (df['mapping_location'].str.extract(pattern))

    return df

def add_gene_info(df):

    pattern = r'(Synt_\d+)'
    df['gene'] = (df['mapping_location'].str.extract(pattern))
    print(df['gene'].unique())

    return df

def pivot_length_table(df):
    df['mapping_location'] = df.index
    df = add_gene_info(df)
    df = add_haplotype_info(df)
    print('all data')
    print(df)
    # Filter out non-unique gene-haplotype combinations
    df = df.groupby(['gene', 'haplotype']).filter(lambda x: len(x) == 1)
    print('filtered data')
    # Display the filtered DataFrame
    print(df)

        
    # Group by 'gene' and 'haplotype', then get the index of the row with the maximum 'ref_length' for each group
    # idx = df.groupby(['gene', 'haplotype'])['ref_length'].idxmax()

    # # Use the indices to filter the DataFrame
    # df = df.loc[idx].reset_index(drop=True)

    df = df.pivot(index='gene', columns='haplotype', values=['ref_length'])
    df.columns = ['_'.join(col) for col in df.columns.values]
    # Drop rows with NaN values
    print(df)
    df =  df[['ref_length_hap1', 'ref_length_hap2', 'ref_length_hap3', 'ref_length_hap4']]
    df.dropna(inplace=True)
    print('NaN values dropped')
    print(df)
    df.fillna(0, inplace=True)
    # Make values in the 'ref_length' columns integers
    df = df.astype(int)
    
    return df


def check_length_values(row, percent):
    min_value = row.min()
    max_value = row.max()

    # Calculate the 1% range
    ten_percent = min_value * percent * 0.01

    # Check if all values are within 1% of each other
    return max_value >= min_value + ten_percent


def add_length_category(df):
    # Select only the relevant columns
    print(df)

    df_subset = df[['ref_length_hap1', 'ref_length_hap2', 'ref_length_hap3', 'ref_length_hap4']]


    df_subset = df_subset.astype(int)
    print(df_subset)

    # print max and min values of each row
    print(df_subset.apply(lambda x: x.max() - x.min(), axis=1))
    df['length_category'] = 'unclassified'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 0), axis=1), 'length_category'] = 'less_1%_difference'
    df.loc[df_subset.apply(lambda x: check_length_values(x, 1), axis=1), 'length_category'] = 'more_1%_difference'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 5), axis=1), 'length_category'] = 'more_5%_difference'

    df.loc[df_subset.apply(lambda x: check_length_values(x, 10), axis=1), 'length_category'] = 'more_10%_difference'


    df.loc[df_subset.apply(lambda x: check_length_values(x, 20), axis=1), 'length_category'] = 'more_20%_difference'

    print(df['length_category'])
    return df

def add_longest_transcript(df):

    df['haplotype_with_longest_annotation'] = df[['ref_length_hap1', 'ref_length_hap2', 'ref_length_hap3', 'ref_length_hap4']].idxmax(axis=1)
    mask = (df[['ref_length_hap1', 'ref_length_hap2', 'ref_length_hap3', 'ref_length_hap4']].nunique(axis=1) == 1)
    df.loc[mask, 'haplotype_with_longest_annotation'] = 'equal_lengths'
    df['haplotype_with_longest_annotation'] = df['haplotype_with_longest_annotation'].str.replace('ref_length_', '')
    #df.loc[df['length_category'] == 'less_1%_difference', 'haplotype_with_longest_annotation'] = 'lengths_within_5%'
    return df

def make_barplot(df):
    plt.figure(figsize=(5, 3))
    # Sort the DataFrame by the 'length_category' column
    custom_order = ['less_1%_difference','more_1%_difference', 'more_5%_difference', 'more_10%_difference', 'more_20%_difference']  # replace with your actual categories
    df['length_category'] = pd.Categorical(df['length_category'], categories=custom_order, ordered=True)
    # chagne the haplotype_with_longest_annotation to a categorical variable
    df['haplotype_with_longest_annotation'] = pd.Categorical(df['haplotype_with_longest_annotation'], categories=['hap2', 'hap4', 'hap1', 'hap3', 'equal_lengths'], ordered=True)
    sns.countplot(x='length_category', hue='haplotype_with_longest_annotation', data=df)
    plt.xlabel('Length Category')
    plt.ylabel('Count')
    # turn the x-axis labels
    plt.xticks(rotation=90)
    plt.title('Counts of Length Categories')
    plt.savefig(f'length_categories.png', dpi=300, bbox_inches='tight', transparent=True)
    plt.close()

    # Also save the data to a file
    df.to_csv(f'length_categories.csv')




def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Get transcript lengths from a GFF file")

    # Add the arguments
    parser.add_argument('GtfFile', metavar='gtf_file', type=str, help='the GFF file to parse')
    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='the GFF file to parse')

    # Parse the arguments
    args = parser.parse_args()

    # read in the gtf file
    df = get_transcript_lengths(args.GtfFile)


    df = pivot_length_table(df)

    transcript_lengths_with_cat = add_length_category(df)
    transcript_lengths_with_cat = add_longest_transcript(transcript_lengths_with_cat)
    # group by lenth category and count
    transcript_lengths_with_cat_group = transcript_lengths_with_cat.groupby(['length_category']).size().reset_index(name='count')
    # save to file 
    #transcript_lengths_with_cat_group.to_csv(f'length_categories.csv', index=False)


    make_barplot(transcript_lengths_with_cat)
    print(transcript_lengths_with_cat)


    

if __name__ == "__main__":
    main()



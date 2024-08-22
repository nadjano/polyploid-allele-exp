import gffutils
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def add_haplotype_info(df, output_prefix):
    # Add a column with the haplotype information
    df['mapping_location'] = df['mapping_location'].astype(str).str.strip()
    if 'RIL' in output_prefix:
        df['haplotype'] = df['mapping_location'].str.split('_').str[0]
    if 'Atlantic' in output_prefix:
        pattern = r'(\dG)'
        df['haplotype'] = df['mapping_location'].str.extract(pattern)
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

def get_transcript_lengths(gff_file, output_prefix):
    # Create a database from the GFF file
    db = gffutils.create_db(gff_file, dbfn='gff.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)

    # Initialize an empty dictionary to store transcript lengths
    transcript_lengths = {}

    # Iterate over all features of type 'transcript' in the GFF file
    for transcript in db.features_of_type('exon'):
        print(transcript)
        if 'Atlantic' in output_prefix:
            if 'Synt' in transcript.id and 'x4' in transcript.id:
                # Calculate the length of the transcript
                length = transcript.end - transcript.start + 1
                
                # Store the length in the dictionary, using the transcript ID as the key
                transcript_lengths[transcript.id] = length
                print(transcript_lengths)
        if 'Orang' in output_prefix:
      
            # Calculate the length of the transcript
            length = transcript.end - transcript.start + 1

            # Store the length in the dictionary, using the transcript ID as the key
            transcript_lengths[transcript.id] = length
    
    df = pd.DataFrame.from_dict(transcript_lengths, orient='index', columns=['ref_length'])
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

def check_length_values(row, percent):
    min_value = row.min()
    max_value = row.max()

    # Calculate the 1% range
    ten_percent = min_value * percent * 0.01

    # Check if all values are within 1% of each other
    return max_value >= min_value + ten_percent


def add_length_category(df, output_prefix):
    # Select only the relevant columns
    print(df)
    if 'Atlantic' in output_prefix:
        df_subset = df[['ref_length_1G', 'ref_length_2G', 'ref_length_3G', 'ref_length_4G']]
    if 'Orang' in output_prefix:
        df_subset = df[['ref_length_hap1', 'ref_length_hap2']]
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


def group_by_gene(df):
    # Group by 'SyntID' and 'haplotype', and sum the 'Length' values
    grouped = df.groupby(['gene', 'haplotype'])['ref_length'].sum()
    
    # Reset the index of the grouped DataFrame
    grouped = grouped.reset_index()
    return(grouped)

def get_haplotype_with_longest_annotation(row):
    if row['ref_length_hap1'] == row['ref_length_hap2']:
        return 'equal_lengths'
    elif row['ref_length_hap1'] > row['ref_length_hap2']:
        return 'ref_length_hap1'
    else:
        return 'ref_length_hap2'    


def add_longest_transcript(df, output_prefix):
    if 'Atlantic' in output_prefix:
        df['haplotype_with_longest_annotation'] = df[['ref_length_1G', 'ref_length_2G', 'ref_length_3G', 'ref_length_4G']].idxmax(axis=1)
        mask = (df[['ref_length_1G', 'ref_length_2G', 'ref_length_3G', 'ref_length_4G']].nunique(axis=1) == 1)
        df.loc[mask, 'haplotype_with_longest_annotation'] = 'equal_lengths'
    if 'Orang' in output_prefix:
        df['haplotype_with_longest_annotation'] = df.apply(get_haplotype_with_longest_annotation, axis=1)
    # Remove the 'ref_length_unique_' prefix
    df['haplotype_with_longest_annotation'] = df['haplotype_with_longest_annotation'].str.replace('ref_length_', '')
    #df.loc[df['length_category'] == 'less_1%_difference', 'haplotype_with_longest_annotation'] = 'lengths_within_5%'
    return df

def pivot_length_table(df, output_prefix):
    df['mapping_location'] = df.index
    df = add_gene_info(df, output_prefix)
    df = add_haplotype_info(df, output_prefix)
    # print genes that appear more than twice in the dataframe
    gene_counts = df['gene'].value_counts()
    genes_more_than_twice = gene_counts[gene_counts > 2]

    print( genes_more_than_twice)
    print(df)
    df = df.pivot(index='gene', columns='haplotype', values=['ref_length'])
    df.columns = ['_'.join(col) for col in df.columns.values]
    # Drop rows with NaN values
    print(df)
    df.dropna(inplace=True)
    df.fillna(0, inplace=True)
    # Make values in the 'ref_length' columns integers
    df = df.astype(int)
    
    return df

def make_barplot(df, output_prefix):
    plt.figure(figsize=(5, 3))
    # Sort the DataFrame by the 'length_category' column
    custom_order = ['less_1%_difference','more_1%_difference', 'more_5%_difference', 'more_10%_difference', 'more_20%_difference']  # replace with your actual categories
    df['length_category'] = pd.Categorical(df['length_category'], categories=custom_order, ordered=True)
    sns.countplot(x='length_category', hue='haplotype_with_longest_annotation', data=df)
    plt.xlabel('Length Category')
    plt.ylabel('Count')
    # turn the x-axis labels
    plt.xticks(rotation=90)
    plt.title('Counts of Length Categories')
    plt.savefig(f'{output_prefix}_length_categories.png', dpi=300, bbox_inches='tight', transparent=True)
    plt.close()


    
def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Get transcript lengths from a GFF file")

    # Add the arguments
    parser.add_argument('GffFile', metavar='gff_file', type=str, help='the GFF file to parse')
    parser.add_argument('output_prefix', metavar='output_prefix', type=str, help='the GFF file to parse')

    # Parse the arguments
    args = parser.parse_args()

    # Get the transcript lengths
    transcript_lengths = get_transcript_lengths(args.GffFile, args.output_prefix)

    print(transcript_lengths)
    # Convert the dictionary to a DataFrame
    
    df = pivot_length_table(transcript_lengths, args.output_prefix)


    print(df)
    transcript_lengths_with_cat = add_length_category(df, args.output_prefix)
    transcript_lengths_with_cat = add_longest_transcript(transcript_lengths_with_cat, args.output_prefix)
    # group by lenth category and count
    transcript_lengths_with_cat_group = transcript_lengths_with_cat.groupby(['length_category']).size().reset_index(name='count')
    # save to file 
    transcript_lengths_with_cat_group.to_csv(f'{args.output_prefix}_length_categories.csv', index=False)


    make_barplot(transcript_lengths_with_cat, args.output_prefix)
    print(transcript_lengths_with_cat)
    

if __name__ == "__main__":
    main()



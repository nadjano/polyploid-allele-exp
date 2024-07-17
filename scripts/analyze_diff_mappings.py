import argparse
import pandas as pd

def add_error_cats(file, output_file):
    df = pd.read_csv(file, sep='\t')
    print(df)
    # Fill na values
    df['mapq_flag_seperate'] = df['mapq_flag_seperate'].fillna('ok')
    df['mapq_flag_competetive'] = df['mapq_flag_competetive'].fillna('ok')
    df['error_category'] = 'no_error'
    df.loc[(df['mapping_category_competetive'].str.contains('unique')) & (df['mapping_category_seperate'].str.contains('unique')) & (df['AS_score_seperate'] > df['AS_score_competetive']), 'error_category'] = 'best_in_comp_not_found'

    df.loc[(df['mapping_category_competetive'].str.contains('unique')) & (df['mapping_category_seperate'].str.contains('unique')) & (df['AS_score_seperate'] < df['AS_score_competetive']), 'error_category'] = 'best_in_sep_not_found'

    df.loc[(df['mapping_category_competetive'].str.contains('unique')) & (df['mapping_category_seperate'].str.contains('multi')), 'error_category'] = 'multi_sep_unique_comp'

    df.loc[(df['mapping_category_competetive'].str.contains('multi')) & (df['mapping_category_seperate'].str.contains('unique')) ,'error_category'] = 'multi_comp_unique_sep'

    df.loc[(df['mapping_category_competetive'].str.contains('multi')) & (df['mapping_category_seperate'].str.contains('multi')) & (df['mapping_location_seperate'].apply(count_commas) > df['mapping_location_competetive'].apply(count_commas)) , 'error_category'] = 'multi_sep_more_comp'

    df.loc[(df['mapping_category_competetive'].str.contains('multi')) & (df['mapping_category_seperate'].str.contains('multi')) & (df['mapping_location_seperate'].apply(count_commas) < df['mapping_location_competetive'].apply(count_commas)) , 'error_category'] = 'multi_comp_more_sep'



    count = df[['error_category']].value_counts()
    print(count)
    # save to file
    df.to_csv(f'{output_file}_error_cats.tsv', sep='\t', index=False)
    print(df)

def count_commas(value):
    return value.count(',') + 1


def group_categories(file, output_file):
    df = pd.read_csv(file, sep='\t')
 
    df_group = df.groupby(['mapping_category_seperate','mapping_category_competetive']).size().reset_index(name='counts')
    print(df)

    df_group_gene = df.groupby(['mapping_location_seperate','mapping_location_competetive']).size().reset_index(name='counts')

    df_group_gene.sort_values(by='counts', ascending=False, inplace=True)
    print(df_group_gene)

    # save to file
    df_group.to_csv(f'{output_file}_grouped_categories.tsv', sep='\t', index=False)

def group_mapping_location(file, output_file):
    df = pd.read_csv(file, sep='\t')
 
    df_group = df.groupby(['mapping_location_seperate','mapq_seperate', 'mapping_location_competetive', 'mapq_competetive']).size().reset_index(name='counts')

    df_group.sort_values(by='counts', ascending=False, inplace=True)
    

    # save to file
    df_group.to_csv(f'{output_file}_grouped_mapping_location.tsv', sep='\t', index=False)

    print(df_group)


def track_reads(file, gene):
    df = pd.read_csv(file, sep='\t')
    #hap1_gene-LOC134761659

    reads_for_gene = df[df['mapping_location_seperate'].str.contains(gene)]

    # save to file
    reads_for_gene.to_csv(f'{gene}_reads.tsv', sep='\t', index=False)


def find_genes(file):
    df = pd.read_csv(file, sep='\t')
    print(df)
    # find the row with the max value for count_seperate
    genes = df.loc[df['count_seperate'].idxmax()]

    # Print the rows with mapping_location == 'hap1_gene-LOC134761659'

    print(df[df['mapping_location'].isin(['hap1_gene-LOC134761659', 'hap2_gene-LOC129052667'])])
   

    return genes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add prefix to fasta sequence IDs.')
    parser.add_argument('-ir', '--input_reads', required=True, help='Input diff file')
    parser.add_argument('-ig', '--input_genes', required=True, help='Input diff file')

    parser.add_argument('-o', '--output', required=True, help='Output prefix')

    args = parser.parse_args()

    
    #group_categories(args.input_reads, args.output)
    add_error_cats(args.input_reads, args.output)
    #group_mapping_location(args.input_reads, args.output)
    #find_genes(args.input_genes)
    #track_reads(args.input_reads, 'hap1_gene-LOC134761659')
    #track_reads(args.input_reads, 'hap2_gene-LOC129052667')
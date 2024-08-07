import pandas as pd
import argparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

pd.options.mode.copy_on_write = True 

def main(args):
    # Read syntelogs file and filter for genes with 4 copies
    syntelogs_df = pd.read_csv(args.syntelogs_file, names=['Syntelog', 'gene_id'], header=None)
    syntleogs_4copies = filter_dataframe_4copies(syntelogs_df)

    # Read GFF3 file
    gff3_df = read_gff3_gzipped_to_df(args.gff3_file)
    
    # Filter syntelogs dataframe to only contain gene_ids that are in the GFF3 file
    filtered_syntelogs, gff_df_mRNA = check_gene_numbers(gff3_df, syntelogs_df)
    print(filtered_syntelogs)
    # Extract haplotypes from gene_ids
    filtered_syntelogs_with_hap = get_haplotyes(filtered_syntelogs)
    print(filtered_syntelogs_with_hap)
    # Filter out syntelogs on S or chrx_0
    filtered_syntelogs_chr_only = remove_unassigned_synt(filtered_syntelogs_with_hap)
    
    # Count haplotypes in syntelogs
    synt_hap_count, synt_group_freq, syntelogs_with_freq = count_haplotypes(filtered_syntelogs_chr_only)
    
    # Rename gene_ids in syntelogs
    rename_gene_ids(syntelogs_with_freq)
    
    # Rename genes in GFF3 file
    #renamed_mRNAgff = rename_genes_in_gff(gff_df_mRNA, syntelogs_with_freq)
    renamed_gff = rename_genes_in_gff(gff3_df, syntelogs_with_freq)
    # Save modified GFF3 file
    save_gff(renamed_gff, args.gff3_outfile)

def read_gff3_gzipped_to_df(file_path):
    """
    Read a gzipped GFF3 file into a pandas DataFrame.
    """
    data = []
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):  # Skip comment lines
                fields = line.strip().split('\t')
                data.append(fields)
 
    # Define column names based on GFF3 format
    columns = ['chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    # Create DataFrame
    df = pd.DataFrame(data, columns=columns)
        
    # Remove Names in attributes
    df['attributes'] = df['attributes'].str.replace(r';Name=Soltu\.Atl_v3.*', '', regex=True)
    # Extract the ID from attribute column
    df['attributes'] = df['attributes'].str.replace('ID=', '',)
    df['attributes'] = df['attributes'].str.replace('Parent=', '',)
    df['attributes'] = df['attributes'].apply(lambda x: remove_last_digit(x))
    print(df)
    return df

def remove_last_digit(s):
    last_dot_index = s.rfind('.')
    return s[:last_dot_index] if last_dot_index != -1 else s

def check_gene_numbers(gff_df, synt_df):
    # Number of unique mRNA ids in GFF3
    gff_df_mRNA = gff_df[gff_df['type'] == "mRNA"]
    gff_df_mRNA['type'] = 'exon'
    
    synt_df['gene_id'] = synt_df['gene_id'].apply(lambda x: remove_last_digit(x))
    print(synt_df)
    # Filter syntelogs dataframe to only contain gene_ids that are in the GFF3 file
    overlapping_syntelogs = synt_df[(synt_df['gene_id'].isin(gff_df_mRNA['attributes']))]
    print("Number of genes not in GFF3 but in syntelogs file:")
    print(len(synt_df['gene_id']) - len(overlapping_syntelogs['gene_id']))
    
    return overlapping_syntelogs, gff_df_mRNA

def filter_dataframe_4copies(df):
    # Group by the 'Synt' column and count the number of unique gene IDs
    synt_counts = df.groupby('Syntelog')['gene_id'].nunique()

    # Filter the DataFrame to only contain rows where the Synt column has four different gene IDs
    filtered_df = df[df['Syntelog'].isin(synt_counts[synt_counts == 4].index)]
    
    return filtered_df

def get_haplotyes(df):
    df["Haplotype"] = df["gene_id"].apply(lambda x: x.split(".")[-1][0:4])
    df['Haplotype'] = df['Haplotype'].replace('^S.*', 'S', regex = True)
    df['Haplotype'] = df['Haplotype'].apply(lambda x: 'chr' + x)
    print(df)
    return df

def rename_gene_ids(df):
    df["Haplotype_Count"] = df["Haplotype_Count"].astype(str)
    df["new_gene_id"] = df["Syntelog"]  + "_" + df["Haplotype"] + "G_x" + df["Haplotype_Count"]
    print(df)
    return df

def rename_genes_in_gff(gff_df, synt_df):
    merged_gff = gff_df.merge(synt_df, how='left', left_on='attributes', right_on='gene_id')

    # Replace 'column1' values with 'column2' values where available
    merged_gff['attributes'] = merged_gff['new_gene_id'].fillna(merged_gff['attributes'])
    merged_gff = merged_gff.drop(columns = ['Syntelog', 'gene_id', 'Haplotype',  'Haplotype_Count', 'new_gene_id'])
    

    # Add "ID" to the 'attributes' column 
    merged_gff['attributes'] = merged_gff['attributes'].apply(lambda x: 'ID=' + x)
    # Add "Parent" to the 'attributes' column if 'type' is 'exon'
    merged_gff.loc[merged_gff['type'] == 'exon', 'attributes'] = merged_gff.loc[merged_gff['type'] == 'exon', 'attributes'].str.replace('ID=', 'Parent=')
    merged_gff.loc[merged_gff['type'] == 'CDS', 'attributes'] = merged_gff.loc[merged_gff['type'] == 'CDS', 'attributes'].str.replace('ID=', 'Parent=')
    
    return merged_gff


def count_haplotypes(df_filtered):
    synt_hap_count = df_filtered.groupby('Syntelog')['gene_id'].nunique().reset_index()
    synt_hap_count.columns = ['Syntelog', 'Haplotype_Count']
    
    synt_group_freq = synt_hap_count.groupby('Haplotype_Count'). \
                                    size().reset_index(name='Frequency')
    synt_group_freq["gene_number"] = synt_group_freq["Haplotype_Count"] * synt_group_freq["Frequency"]
    synt_hap_count = pd.DataFrame(synt_hap_count)

    synt_hap_count['Syntelog'] = synt_hap_count['Syntelog'].astype(str)
    df_filtered['Syntelog'] = df_filtered['Syntelog'].astype(str)

    df_joined = pd.merge(df_filtered, synt_hap_count, on='Syntelog')
    print(df_joined)
    return synt_hap_count, synt_group_freq, df_joined

def remove_unassigned_synt(synt_df):
    filtered_df = synt_df[~synt_df['Haplotype'].str.contains('_0')]
    filtered_df = filtered_df[~filtered_df['Haplotype'].str.contains('S')]
    return filtered_df

def save_gff(gff_df, output_gff):
    gff_df.to_csv(output_gff, sep='\t', header=False, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Description of your script")
    
    # Add arguments here using add_argument() method
    parser.add_argument("-s", "--syntelogs_file", 
                        dest="syntelogs_file",
                        required=True,
                        help="File containing syntelog groups")

    parser.add_argument("-g", "--gff3", 
                            dest="gff3_file", 
                            required=True, 
                            help="File containing gene annotation of all genes")
                            #
    parser.add_argument("-o", "--gff3_out", 
                            dest="gff3_outfile", 
                            required=True, 
                            help="Adapted gff3 file with only mRNA fields and updated gene_id names for syntelogs")
    # Add more arguments as needed
    
    args = parser.parse_args()
    
    main(args)

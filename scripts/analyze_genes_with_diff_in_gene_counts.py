
import pandas as pd
import argparse

def main(file_path):
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t')

    # Display the DataFrame
 
    # Filter to only include genes with the specified ratio conditions
    df_filtered = df[(df['count_seperate'] / (df['count_seperate'] + df['count_competetive']) > 0.7) | 
                 (df['count_seperate'] / (df['count_seperate'] + df['count_competetive']) < 0.3)]
    #df_filtered = df_filtered[(df_filtered['count_seperate'] + df_filtered['count_competetive']) > 20]

    # Extract gene name from 'mapping_location' and create a new column 'gene_name'
    df_filtered['gene_name'] = df_filtered['mapping_location'].apply(lambda x: x.split('-')[1])

    # Sort the DataFrame by 'gene_name'
    df_sorted = df_filtered.sort_values(by='gene_name')
    # Filter to only incldue genes that have more than 20 counts per gene
    df_sorted = df_sorted[(df_sorted['count_seperate'] + df_sorted['count_competetive']) > 40]
    # save the sorted DataFrame to a new TSV file
    #df_sorted = df_sorted[df_sorted.duplicated(subset=['gene_name'], keep=False)]
    df_sorted.to_csv('S4_Orangutan_P_gene_diff_sorted.tsv', sep='\t', index=False)
    # print number of unique gene names
    print(df_sorted['gene_name'].nunique())

    # Display the sorted DataFrame
    print(df_sorted)
    # print where gene name is not unique
    

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Read a TSV file into a pandas DataFrame.')
    parser.add_argument('file_path', type=str, help='Path to the TSV file')

    # Parse the arguments
    args = parser.parse_args()

    # Call the main function with the provided file path
    main(args.file_path)



# python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/analyze_genes_with_diff_in_gene_counts.py /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/out_ms/scatter_plots/Orangutan_P_tsv/Orangutan_P_gene_diff.tsv
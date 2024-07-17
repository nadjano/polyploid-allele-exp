import pandas as pd
import argparse

def read_file(file_path):
    # Read file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t')

    # Print the DataFrame
    print(df)

    # Count the number of unique IDs
    unique_ids = df['qseqid'].nunique()
    print(f'Number of unique IDs: {unique_ids}')

    # Count the number of IDs that have exact duplications
    duplicated_ids = df[df.duplicated(subset=['qseqid'], keep=False)].shape[0]
    print(f'Number of duplicated IDs: {duplicated_ids}')
    return df

def get_primary_mapping(df):
    df['bitscore'] = df['bitscore'].astype(float)
    # If there are multiple mappings per read with the same AS score, drop them
    # Get the 'read_id's to remove
    df_max_filter, num_multi_max = filter_max_score(df)
    
    print("Number of multimappers", num_multi_max)
    print("number of mappings", len(df_max_filter))
    return df_max_filter


def filter_max_score(df):
    # Get rows where 'AS_score' is maximum for each 'read_id'
    idxmax = df.groupby('qseqid')['bitscore'].transform('max') == df['bitscore']
    df_max = df.loc[idxmax]

    # Get rows where 'read_id' and 'AS_score' are duplicated
    non_unique_max = df_max.duplicated(subset=['qseqid', 'bitscore'], keep=False)
    df_multi = df_max[non_unique_max]

    # Get the unique read ids in the multi-mapped reads
    df_multi_unique = df_multi['qseqid'].unique()

    # Filter out the non-unique max scores
    df_max_filter = df_max[~non_unique_max]

    # Return the filtered DataFrame and the number of reads with multiple max scores
    return df_max_filter, len(df_multi_unique)

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Read a tab-separated file into a DataFrame")

    # Add an argument for the file path
    parser.add_argument('file_path', help="The path to the file to read")

    # Parse the arguments
    args = parser.parse_args()

    # Read the file
    df = read_file(args.file_path)
    df_primiary = get_primary_mapping(df)

if __name__ == "__main__":
    main()
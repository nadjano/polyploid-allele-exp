import argparse
import pandas as pd

pd.options.mode.chained_assignment = None

def parse_line(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
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
    return read_id, mapping_location, mapq, AS_score, NM_score

def paf_to_dataframe(file):
    data = [parse_line(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id', 'mapping_location', 'mapq', 'AS_score', 'NM_score'])
    # return a test df
    #df = df.head(1000)
    return df

def filter_max_score_keep_dups(df):
    # Fill na values with 0
    df['AS_score'] = df['AS_score'].fillna(0)
    # Get rows where 'AS_score' is maximum for each 'read_id' 
    idxmax = df.groupby('read_id')['AS_score'].transform('max') == df['AS_score']

    df_max = df.loc[idxmax]

    # where 'read_id', 'mapping location' and 'AS_score' are duplicated select one of them
    non_unique_max = df_max.duplicated(subset=['read_id', 'mapping_location', 'AS_score'])

    # Remove all duplicated rows
    # df_max = df_max.drop_duplicates(subset=['read_id', 'mapping_location', 'AS_score'], keep="last")
    print(df_max)
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

def get_primary_mapping(df):
    df = df[df['mapping_location'] != '*']
    df['mapq'] = df['mapq'].astype(int)
    df = df[df['mapq'] >= 1]
    df['AS_score'] = df['AS_score'].astype(float)
    # If there are multiple mappings per read with the same AS score, drop them
    # Get the 'read_id's to remove
    df_max_filter, num_multi_max = filter_max_score(df)

    return df_max_filter

def compare_mappings(df1, df2, score):
    merged = pd.merge(df1, df2, on='read_id', suffixes=('_separate', '_competetive'))
    print(merged['mapping_location_separate'])
    print(merged['mapping_location_competetive'])
    same_mapping = (merged['mapping_location_separate'] == merged['mapping_location_competetive']).sum()
    different_mapping = (merged['mapping_location_separate'] != merged['mapping_location_competetive']).sum()
    # Save the rows that have different mapping locations
    print(merged[merged['mapping_location_separate'] != merged['mapping_location_competetive']])

    same_score = ((merged['mapping_location_separate'] != merged['mapping_location_competetive']) &
                  (merged[score + '_separate'] == merged[score + '_competetive'])).sum()
    different_score = ((merged['mapping_location_separate'] != merged['mapping_location_competetive']) &
                       (merged[score + '_separate'] != merged[score + '_competetive'])).sum()
    return same_mapping, different_mapping, same_score, different_score

def pivot_table(df):
    df = df.sort_values(by=['read_id', 'mapping_location'])
    return df.pivot_table(index=['read_id'], values=['mapping_location', 'mapq', 'AS_score'], aggfunc=lambda x: ','.join(str(v) for v in x)) # ['sum', 'count']

def check_same_values(value):
    parts = value.split(',')
    return len(set(parts)) == 1

def join_dfs(dataframes):
    for i in range(len(dataframes)):
        dataframes[i] = dataframes[i].sort_values(by=['read_id', 'mapping_location'])
        dataframes[i] = dataframes[i][['read_id', 'mapping_location', 'AS_score', 'NM_score', 'mapq']]
        dataframes[i] = filter_max_score_keep_dups(dataframes[i])
        dataframes[i] = pivot_table(dataframes[i])

    merged = pd.merge(dataframes[0], dataframes[1], on='read_id', suffixes=('_hap1', '_hap2'), how= 'outer')
    merged = pd.merge(merged, dataframes[2], on='read_id', suffixes=('_minimap', 'competetive'), how= 'outer')
    print(merged)
    return merged




def main( paf_comp):
    df = pd.read_csv(paf_comp, sep='\t')
    # Filter for reads that map to Syntelogs
    df = df[df['mapping_location'].str.contains('Synt')]
    df_max = filter_max_score_keep_dups(df)
    pivot_df = pivot_table(df_max)
    print(pivot_df)
    # get frequencies of mapping locations
    mapping_location_freq = pivot_df['mapping_location'].value_counts()
    print(mapping_location_freq)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process one PAF file.')
    parser.add_argument('-i','--paf_competetive', type=str, help='Path to the second TSV file')
    #parser.add_argument('--summary', type=str, help='Path to the summary file')
    #parser.add_argument('--stats', type=str, help='Path to the stats file')
    args = parser.parse_args()
    main(args.paf_competetive)
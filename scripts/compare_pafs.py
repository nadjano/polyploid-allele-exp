import argparse
import pandas as pd

pd.options.mode.chained_assignment = None

def parse_line(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
    read_length = int(fields[1])
    if '_' in fields[5]:
        mapping_location = fields[5].split('_')[1]
    else:
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
    return read_id, read_length, mapping_location, mapq, AS_score, NM_score

def paf_to_dataframe(file):
    data = [parse_line(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id','read_length', 'mapping_location', 'mapq', 'AS_score', 'NM_score'])
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

def filter_max_score(df):
    # Get rows where 'AS_score' is maximum for each 'read_id'
    idxmax = df.groupby('read_id')['AS_score'].transform('max') == df['AS_score']
    df_max = df.loc[idxmax]

    # Get rows where 'read_id' and 'AS_score' are duplicated
    non_unique_max = df_max.duplicated(subset=['read_id', 'AS_score'], keep=False)
    df_multi = df_max[non_unique_max]

    # Get the unique read ids in the multi-mapped reads
    df_multi_unique = df_multi['read_id'].unique()

    # Filter out the non-unique max scores
    df_max_filter = df_max[~non_unique_max]

    # Return the filtered DataFrame and the number of reads with multiple max scores
    return df_max_filter, len(df_multi_unique)

def get_mapping_stats(df1, df2):
    stats = []
    for df in [df1, df2]:
        # Get the number of reads mapped with mapq >5
        df = df[(df['mapping_location'] != '*') & (df['mapq'].astype(int) >=1)]
        total_mappings = len(df)
        reads_mapped = df['read_id'].nunique()
        # Get the number of reads that have the same mapping location and highest AS score
        df_max_filter, num_multi_max = filter_max_score(df)
        print(num_multi_max)
        filtred_unique_mappings = len(df_max_filter)
        stats.append([total_mappings, reads_mapped, num_multi_max, filtred_unique_mappings])
    return pd.DataFrame(stats, columns=['total_mappings', 'unique_mappings', 'dup_primary_mappings', 'filtered_unique_mappings'],  index=['separate', 'competetive'])

def pivot_table(df):
    return df.pivot_table(index=['read_id'], values=['mapping_location','read_length', 'mapq', 'AS_score'], aggfunc=lambda x: ','.join(str(v) for v in x)) # ['sum', 'count']

def check_same_values(value):
    # Return wrong if the value is not contained
    if ',' not in value:
        return False
    else:
        parts = value.split(',')
        return len(set(parts)) == 1

def join_dfs(dataframes):
    for i in range(len(dataframes)):
        dataframes[i] = dataframes[i].sort_values(by=['read_id', 'mapping_location'])
        dataframes[i] = dataframes[i][['read_id','read_length', 'mapping_location', 'AS_score', 'NM_score', 'mapq']]
        dataframes[i] = filter_max_score_keep_dups(dataframes[i])
        dataframes[i] = pivot_table(dataframes[i])

    merged = pd.merge(dataframes[0], dataframes[1], on='read_id', suffixes=('_hap1', '_hap2'), how= 'outer')
    merged = pd.merge(merged, dataframes[2], on='read_id', suffixes=('_minimap', 'competetive'), how= 'outer')
    print(merged)
    return merged

def remove_duplicates(s):
    return ','.join(set(s.split(',')))



def add_mapping_category(merged_df, summary_file, stats_file):
    # set all the mapping_category tags to 'uncalssified'
    merged_df['mapping_category'] = 'unclassified'

    # get the multi-mapped reads that have the same mapping location
    merged_df.loc[(merged_df['AS_score'].apply(check_same_values)) & (~merged_df['mapq'].str.contains("0,0").fillna(False)), 'mapping_category'] = 'multi_comp_error'
    merged_df.loc[(merged_df['mapping_location_hap2'] == merged_df['mapping_location_hap1']) & 
              (merged_df['AS_score_hap2'] == merged_df['AS_score_hap1']) & (merged_df['mapping_location'] == '*') | (merged_df['mapping_location'].isna()) & 
              (~merged_df['mapq'].str.contains("0,0").fillna(False)), 
              'mapping_category'] = 'strange_multi_sep_unique_comp'
    merged_df.loc[(merged_df['mapping_location_hap1'] ==
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap1'] == merged_df['AS_score_hap2']) & (merged_df['mapq'] == "0,0,0,0"), 'mapping_category'] = 'multi_same_2_genes'
    merged_df.loc[(merged_df['mapping_location_hap1'] ==
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap1'] == merged_df['AS_score_hap2']) & (merged_df['mapping_location'] == '*') | (merged_df['mapping_location'].isna()), 'mapping_category'] = 'multi_same_gene_not_comp'
    merged_df.loc[(merged_df['mapping_location_hap1'] ==
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap1'] == merged_df['AS_score_hap2']) & (merged_df['mapping_location'].apply(check_same_values)) & (merged_df['mapq'] == "0,0"), 'mapping_category'] = 'multi_same_1_gene'
    merged_df.loc[(merged_df['mapping_location_hap1'] == merged_df['mapping_location_hap2']) & 
              (merged_df['AS_score_hap1'] == merged_df['AS_score_hap2']) & 
              (merged_df['mapq'].str.contains("0,0,0,0,0")), 'mapping_category'] = 'multi_same_more_than_2_genes'
    merged_df.loc[(merged_df['mapping_location_hap1'] ==
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap1'] > merged_df['AS_score_hap2']) & (merged_df['mapping_location'] == (merged_df['mapping_location_hap1'])) & (merged_df['mapq'] != "0,0"),  'mapping_category'] = 'unique_hap1_better_same_gene'
    merged_df.loc[(merged_df['mapping_location_hap1'] ==
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap2'] > merged_df['AS_score_hap1']) & (merged_df['mapping_location'] == (merged_df['mapping_location_hap2'])) & (merged_df['mapq'] != "0,0"),  'mapping_category'] = 'unique_hap2_better_same_gene'
    merged_df.loc[(merged_df['mapq'].isna()) & (merged_df['AS_score_hap2'] > merged_df['AS_score_hap1']) & (merged_df['mapping_location_hap2'] == (merged_df['mapping_location_hap1'])),  'mapping_category'] = 'unique_hap2_better_not_competetive'
    merged_df.loc[(merged_df['mapq'].isna()) & (merged_df['AS_score_hap1'] > merged_df['AS_score_hap2']) & (merged_df['mapping_location_hap2'] == merged_df['mapping_location_hap1']),  'mapping_category'] = 'unique_hap1_better_not_competetive'
    merged_df.loc[(merged_df['mapq'].isna()) & (merged_df['AS_score_hap1'] == merged_df['AS_score_hap2']) & (merged_df['mapping_location_hap2'] == merged_df['mapping_location_hap1']),  'mapping_category'] = 'multi_not_competetive'
    merged_df.loc[(merged_df['AS_score_hap1'] == merged_df['AS_score']) & (merged_df['mapping_location_hap2'].isna()) & (merged_df['mapping_location'] == merged_df['mapping_location_hap1']), 'mapping_category'] = 'unique_hap1_not_hap2_competetive'
    merged_df.loc[(merged_df['AS_score_hap2'] == merged_df['AS_score']) & (merged_df['mapping_location_hap1'].isna()) & (merged_df['mapping_location'] == merged_df['mapping_location_hap2']), 'mapping_category'] = 'unique_hap2_not_hap1_competetive'
    merged_df.loc[(merged_df['mapping_location_hap1'] !=
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap1'] > merged_df['AS_score_hap2']) & (merged_df['mapping_location'] == (merged_df['mapping_location_hap1'])) & (merged_df['mapq'] != "0,0") & (merged_df['mapping_location_hap2'] != "*"),  'mapping_category'] = 'unique_hap1_better_diff_gene'
    merged_df.loc[(merged_df['mapping_location_hap1'] !=
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap2'] > merged_df['AS_score_hap1']) & (merged_df['mapping_location'] == (merged_df['mapping_location_hap2'])) & (merged_df['mapq'] != "0,0") & (merged_df['mapping_location_hap1'] != "*"),  'mapping_category'] = 'unique_hap2_better_diff_gene'
    merged_df.loc[(merged_df['mapping_location_hap1'] ==
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap1'] > merged_df['AS_score_hap2'])& (merged_df['mapq_hap2'].str.contains("0,0")) & (merged_df['mapping_location'] == (merged_df['mapping_location_hap1'])) & (merged_df['mapq'].str.contains("0,0")),  'mapping_category'] = 'multi_hap1_better_same_genes'
    merged_df.loc[(merged_df['mapping_location_hap1'] ==
    merged_df['mapping_location_hap2']) & (merged_df['AS_score_hap2'] > merged_df['AS_score_hap1'])& (merged_df['mapq_hap2'].str.contains("0,0")) & (merged_df['mapping_location'] == (merged_df['mapping_location_hap1'])) & (merged_df['mapq'].str.contains("0,0")),  'mapping_category'] = 'multi_hap2_better_same_genes'
    merged_df.loc[(merged_df['mapping_location'] == merged_df['mapping_location_hap1']) & (merged_df['mapping_location_hap2'] == "*") & (merged_df['mapping_location'] != "*") & (merged_df['mapq']!= "0,0"), 'mapping_category'] = 'unique_hap1_not_hap2'
    merged_df.loc[(merged_df['mapping_location'] == merged_df['mapping_location_hap2']) & (merged_df['mapping_location_hap1'] == "*") & (merged_df['mapping_location'] != "*") & (merged_df['mapq']!= "0,0"), 'mapping_category'] = 'unique_hap2_not_hap1'
    merged_df.loc[(merged_df['mapping_location'] == merged_df['mapping_location_hap2']) & (merged_df['mapping_location_hap1'] == merged_df['mapping_location_hap2']) & (merged_df['mapq'] == "0"), 'mapping_category'] = 'multi_comp_strange'


    merged_df.loc[((merged_df['mapping_location_hap1'] == '*') | (merged_df['mapping_location_hap1'].isna())) & ((merged_df['mapping_location_hap2'] == '*') | (merged_df['mapping_location_hap2'].isna())) & ((merged_df['mapping_location'] == '*') | (merged_df['mapping_location'].isna())), 'mapping_category'] = 'unmapped'

    merged_df =merged_df.applymap(remove_duplicates)
    merged_df['read_length'] = merged_df['read_length'].astype(int)

    # Filter for min_read length
    merged_df_lt_200 = merged_df[merged_df['read_length'] >= 200]
    # save the merged df
    merged_df.to_csv(summary_file, sep = '\t')
    # get the mapping category stats
    mapping_category_stats = merged_df['mapping_category'].value_counts()
    mapping_category_lt_200_stats = merged_df_lt_200['mapping_category'].value_counts()
    print(mapping_category_stats)
    mapping_category_stats.to_csv(stats_file, sep = '\t')
    mapping_category_lt_200_stats.to_csv(f'{stats_file}_lt_200', sep = '\t')
    return merged_df

def main(paf1, paf2, paf_comp, summary, stats_file):
    dfs = [paf_to_dataframe(paf) for paf in [paf1, paf2, paf_comp]]
    for df in dfs:
        print("Unique read ids" + str(df['read_id'].nunique()))
    merged = join_dfs(dfs)
    with_cat = add_mapping_category(merged, summary, stats_file)
    print(with_cat)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process two PAF files.')
    parser.add_argument('--paf_seperate1', type=str, help='Path to the first seperate PAF file')
    parser.add_argument('--paf_seperate2', type=str, help='Path to the second separate PAF file')
    parser.add_argument('--paf_competetive', type=str, help='Path to the second PAF file')
    parser.add_argument('--summary', type=str, help='Path to the summary file')
    parser.add_argument('--stats', type=str, help='Path to the stats file')
    args = parser.parse_args()
    main(args.paf_seperate1, args.paf_seperate2, args.paf_competetive, args.summary, args.stats)
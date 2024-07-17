import argparse
import pandas as pd
import re

# For testing: 
### python /blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison/scripts/get_multimappers_from_paf.py --experiment Atlantic --paf_seperate  test_data/align_seperate_SRR14993893_test.paf  --paf_competetive test_data/align_competetive_SRR14993893_test.paf --blast test_data/Atlantic_SRR14993893_head.BLAST_test.tsv  --summary summ.test --mapping mapping.test

pd.options.mode.chained_assignment = None

def parse_line_paf(line):
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
    # Only select first 500 lines for testing
    #data = [parse_line_paf(line) for line in open(file, 'r')]
    # Only select first 500 lines for testing
    with open(file, 'r') as f:
        data = [parse_line_paf(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id', 'mapping_location', 'mapq', 'AS_score', 'NM_score'])
    #df = df[1:1000]
    return df

def blast_to_dataframe(file):
    df = pd.read_csv(file, sep='\t', header=0)
    # remove the "@" from the start of the qseqid
    df['qseqid'] = df['qseqid'].str[1:]
    # replace qseq_id with read_id
    df = df.rename(columns={'qseqid': 'read_id', 'mismatch': 'NM_score', 'bitscore': 'AS_score', 'stitle': 'mapping_location'})

    # Get a subset for testing
    #df = df[1:1000]
    return df


def add_haplotype_info(df, experiment):
    # Add a column with the haplotype information
    df['mapping_location'] = df['mapping_location'].astype(str).str.strip()
    if experiment == 'RIL':
        df['haplotype'] = df['mapping_location'].str.split('_').str[0]
    if experiment == 'RIL_updated':
        df['haplotype'] = df['mapping_location'].str.split('_').str[0]
    elif experiment == 'Atlantic':
        pattern = r'(\dG)'
        df['haplotype'] = df['mapping_location'].str.extract(pattern)
        print(df['haplotype'])

    return df

def get_mapping_stats(df1, df2, df3):
    stats = []
    for df in [df1, df2, df3]:
        # Get the number of reads mapped with mapq >5
        df = df[(df['mapping_location'] != '*')] #& (df['mapq'].astype(int) >=1)]
        total_mappings = len(df)
        reads_mapped = df['read_id'].nunique()
        # Get the number of reads that have the same mapping location and highest AS score
        df_max_filter, num_multi_max = filter_max_score(df)
        print(num_multi_max)
        filtred_unique_mappings = len(df_max_filter)
        stats.append([total_mappings, reads_mapped, num_multi_max, filtred_unique_mappings])
    stats_df = pd.DataFrame(stats, columns=['total_mappings', 'unique_read_ids', 'dup_primary_mappings', 'filtered_unique_mappings'],  index=['mm2_separate', 'mm2_competetive', 'blast'])
    return stats_df

def save_stats_to_tsv(df, outfile):
    # Calculate the total number of reads
    total_reads = len(df)
    total_reads_all = (df[~df['mapping_location'].isna() & ~df['mapping_location_separate'].isna() & ~df['mapping_location_competetive'].isna()])
    total_reads_all_len = len(total_reads_all)
    df_only_single = total_reads_all[total_reads_all['multi_mapping_type'] == 'single_best']
    df_only_multi = total_reads_all[total_reads_all['multi_mapping_type'].str.contains('multi')]
    df_only_multi_diff = total_reads_all[total_reads_all['multi_mapping_type'] == 'multi_different_haplotype']
    # Calculate the statistics
    stats = {
        'total_reads': total_reads,
        'reads_mapped_for_all': total_reads_all_len,
        'reads_mapped_mm_separate': len(df[~df['mapping_location_separate'].isna()]),
        'reads_mapped_mm_competetive': len(df[~df['mapping_location_competetive'].isna()]),
        'reads_mapped_blast': len(df[~df['mapping_location'].isna()]),
        'reads_mapped_mm_separate_type': df['multi_mapping_type_separate'].value_counts().to_dict(),
        'reads_mapped_mm_competetive_type': df['multi_mapping_type_competetive'].value_counts().to_dict(),
        'reads_mapped_blast_type': df['multi_mapping_type'].value_counts().to_dict(),
        'agreement_blast_intersect': total_reads_all['mapping_location_agreement_blast'].sum(),
        'agreement_mm_intersect': total_reads_all['mapping_location_agreement_mm'].sum(),
        'agreement_all_intersect': total_reads_all['mapping_location_agreement_all'].sum(),
        # Get agreements for the different categories
        'agreement_all_percentage_intersect_only_single_best': df['mapping_location_agreement_blast'].sum() / len(df) * 100,
        'agreement_all_percentage_intersect_only_single_best': df_only_single['mapping_location_agreement_blast'].sum() / len(df_only_single) * 100,
        'agreement_all_percentage_intersect_multi_hap_best': df_only_multi['mapping_location_agreement_blast'].sum() / len( df_only_multi) * 100    
    }

    # Convert the dictionary to a DataFrame
    stats_df = pd.DataFrame.from_dict(stats, orient='index')

    # Save the DataFrame to a TSV file
    stats_df.to_csv(outfile, sep='\t', header=False)


def calculate_multimapping_frequencies(df):
    # Calculate the frequencies of the multi-mapping types
    frequencies = df['multi_mapping_type'].value_counts()
    return frequencies

def get_primary_mapping(df):

    df = df[df['mapping_location'] != '*']

    df['AS_score'] = df['AS_score'].astype(float)
    
    df_max_filter = filter_max_score_keep_dups(df)

    return df_max_filter

def filter_max_score_keep_dups(df):
    # Get rows where 'AS_score' is maximum for each 'read_id'
    idxmax = df.groupby('read_id')['AS_score'].transform('max') == df['AS_score']

    df_max = df.loc[idxmax]

    # where 'read_id', 'mapping location' and 'AS_score' are duplicated select one of them
    non_unique_max = df_max.duplicated(subset=['read_id', 'mapping_location', 'AS_score'])
    # For the first duplicated row set add a tag 
    df_max.loc[non_unique_max, 'duplicated_in_gene'] = "yes"
    df_max.loc[~non_unique_max, 'duplicated_in_gene'] = "no"
    # Remove all duplicated rows
    df_max = df_max.drop_duplicates(subset=['read_id', 'mapping_location', 'AS_score'], keep="last")
    print(df_max)
    # Return the filtered DataFrame and the number of reads with multiple max scores
    return df_max

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

def add_multi_mapping_type_tag(df):
    # for the pivoted dataframe add a tag for the type of multi mapping
    df['multi_mapping_type'] = None
    # if the mapping location contains no "," then it is a single best mapping"
    df.loc[~df['mapping_location'].str.contains(','), 'multi_mapping_type'] = 'single_best'
    # If the mapping location contains ',' then check if the haplotypes are the same (separated by ',')
    df.loc[df['mapping_location'].str.contains(','), 'multi_mapping_type'] = df['haplotype'].apply(lambda x: 'multi_same_haplotype' if len(set(x.split(','))) == 1 else 'multi_different_haplotype')
    percentage = calculate_multimapping_frequencies(df)
    print(percentage)
    return df

def joindfs(df1, df2, df3, experiment):
    # Sort tables on read_id and mapping location

    dataframes = [df1, df2, df3]

    for i in range(len(dataframes)):
        dataframes[i] = dataframes[i].sort_values(by=['read_id', 'mapping_location'])
        dataframes[i] = dataframes[i][['read_id', 'mapping_location', 'AS_score', 'NM_score', 'duplicated_in_gene']]
        dataframes[i] = add_haplotype_info(dataframes[i], experiment)
        dataframes[i] = pivot_table(dataframes[i])
        dataframes[i] = add_multi_mapping_type_tag(dataframes[i])
        print(dataframes[i])

    merged = pd.merge(dataframes[0], dataframes[1], on='read_id', suffixes=('_separate', '_competetive'), how= 'inner')
    merged = pd.merge(merged, dataframes[2], on='read_id', suffixes=('_minimap', '_blast'),  how= 'inner')
    print(merged)
    return merged

def pivot_table(df):
    return df.pivot_table(index=['read_id'], values=['mapping_location', 'duplicated_in_gene', 'haplotype'], aggfunc=lambda x: ','.join(x) ) # ['sum', 'count']

def add_mapping_location_agreement_tag(df):
    # set a tag when all mapping locations are the same

    df['mapping_location_agreement_all'] = (df['mapping_location_separate'] == df['mapping_location_competetive']) & (df['mapping_location_competetive'] == df['mapping_location'])
    df['mapping_location_agreement_mm'] = (df['mapping_location_separate'] == df['mapping_location_competetive'])
    df['mapping_location_agreement_blast'] = (df['mapping_location_competetive'] == df['mapping_location'])


    return df

def get_percentages(df, results_file, name):
    # Get the percentage of reads that have mapping_location_agreement_blast, mapping_location_agreement_mm, mapping_location_agreement
    total_reads = len(df)
    # Get the percentage of reads that have mapping_location_agreement_blast, mapping_location_agreement_mm, mapping_location_agreement
    # Avid dividing by zero
    if total_reads == 0:
        return 0, 0, 0
    agreement_blast = df['mapping_location_agreement_blast'].sum() / total_reads
    agreement_mm = df['mapping_location_agreement_mm'].sum() / total_reads
    agreement = df['mapping_location_agreement_all'].sum() / total_reads

    # Write the results to a file
    with open(results_file, 'a') as f:
        f.write(f'{name}\n')
        f.write(f'Total reads: {total_reads}\n')
        f.write(f'Agreement between blast and minimap: {agreement_blast}\n')
        f.write(f'Agreement between minimap and minimap: {agreement_mm}\n')
        f.write(f'Agreement between all three methods: {agreement}\n')

    return agreement_blast, agreement_mm, agreement


def main(paf1, paf2, blast, summary, mapping, experiment):
    dfs = [paf_to_dataframe(paf) for paf in [paf1, paf2]]
    df_blast = blast_to_dataframe(blast)

    # Get the mapping stats for the three different mapping methods
    mapping_stats = get_mapping_stats(dfs[0], dfs[1], df_blast)
    # Save the stats to a file
    mapping_stats.to_csv(summary + '_stats', sep='\t')

    # Filter for max AS score
    primary_mappings_mm = [get_primary_mapping(df) for df in dfs]
    primary_mappings_blast = get_primary_mapping(df_blast)

    # Merge the dataframes
    merged = joindfs(primary_mappings_mm[0], primary_mappings_mm[1], primary_mappings_blast, experiment)

    with_tag = add_mapping_location_agreement_tag(merged)
    # Call the function
    save_stats_to_tsv(with_tag, summary)
    # Save the data to a file
    with_tag.to_csv(mapping, sep='\t')
    # Get the percentages of reads that have the same mapping location
    # Filter to only include that have a , in all of the mapping locations --> more than one best mapping location in all methods
    with_tag_all_dup = with_tag[with_tag['mapping_location_separate'].str.contains(',') & with_tag['mapping_location_competetive'].str.contains(',') & with_tag['mapping_location'].str.contains(',')]
    # Filter to only include that have at least one , in one of the mapping locations --> more than one best mapping location
    with_tag_only_dup = with_tag[with_tag['mapping_location_separate'].str.contains(',') | with_tag['mapping_location_competetive'].str.contains(',') | with_tag['mapping_location'].str.contains(',')]
    # Filter to only include that no ',' in any of the mapping locations --> only one best mapping location
    with_tag_only_single_best = with_tag[~(with_tag['mapping_location_separate'].str.contains(',') | with_tag['mapping_location_competetive'].str.contains(',') | with_tag['mapping_location'].str.contains(','))]
    # Save the data to a file
    summary = summary + '_agreement'
    get_percentages(with_tag, summary, 'all_intersect_mapped_reads')
    get_percentages(with_tag_only_dup, summary, '_some_dupplicated_best_match')
    get_percentages(with_tag_only_single_best, summary, '_only_single_best')
    get_percentages(with_tag_all_dup, summary, '_all_dupplicated_best_match')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process two PAF files.')
    parser.add_argument('--experiment', type=str, help='Name of the experiment')
    parser.add_argument('--paf_seperate', type=str, help='Path to the first PAF file')
    parser.add_argument('--paf_competetive', type=str, help='Path to the second PAF file')
    parser.add_argument('--blast', type=str, help='Path to the blast file')
    parser.add_argument('--summary', type=str, help='Path to the summary file')
    parser.add_argument('--mapping', type=str, help='Path to the mapping file')
    # parser.add_argument('--stats', type=str, help='Path to the stats file')
    args = parser.parse_args()
    main(args.paf_seperate, args.paf_competetive, args.blast, args.summary, args.mapping, args.experiment)
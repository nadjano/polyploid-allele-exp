import argparse
import pandas as pd

pd.options.mode.chained_assignment = None

def parse_line(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
    read_length = int(fields[1])
    mapping_location = fields[5]
    mapq = fields[11]
    AS_score = None
    NM_score = None
    for field in fields[12:]:
        if field.startswith('NM'):
            NM_score = field.split(':')[2]
        if field.startswith('ms'):
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
    # Convert 'AS_score' to integers
    df['AS_score'] = df['AS_score'].astype(int)
    # Convert 'mapq' to integers
    df['mapq'] = df['mapq'].astype(int)
    # Get rows where 'AS_score' is maximum for each 'read_id' 
    idxmax = df.groupby('read_id')['AS_score'].transform('max') == df['AS_score']

    df_max = df.loc[idxmax]

    # where 'read_id', 'mapping location' and 'AS_score' are duplicated select one of them
    non_unique_max = df_max.duplicated(subset=['read_id', 'mapping_location', 'AS_score'])

    # Remove all duplicated rows
    # df_max = df_max.drop_duplicates(subset=['read_id', 'mapping_location', 'AS_score'], keep="last")
    print(df_max)
    # Add haplotype information
    # Group by read_id and sum other columns
    df_max =  df_max.pivot_table(index=['read_id'], values=['mapping_location','mapq', 'AS_score', 'read_length', 'haplotype'], aggfunc=lambda x: ','.join(str(v) for v in x))
    print(df_max[['mapping_location']])

    # Remove duplicates in read_length
    df_max['read_length'] = pd.to_numeric(df_max['read_length'].apply(remove_duplicates), errors='coerce')
    df_max['mapq'] = pd.to_numeric(df_max['mapq'].apply(remove_duplicates), errors='coerce')
    df_max = df_max.dropna(subset=['read_length'])
    df_max['read_length'] = df_max['read_length'].astype(int)
    # Fill na values with 0
    df_max['mapq'] = df_max['mapq'].fillna(0)
    df_max['mapq'] = df_max['mapq'].astype(int)
    return df_max


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



def pivot_table(df):
    return df.pivot_table(index=['read_id'], values=['mapping_location','read_length', 'mapq', 'AS_score'], aggfunc=lambda x: ','.join(str(v) for v in x)) # ['sum', 'count']

def check_same_values(value):
    # Return wrong if the value is not contained
    if ',' not in value:
        return False
    else:
        parts = value.split(',')
        return len(set(parts)) == 1


def remove_duplicates(s):
    return ','.join(set(s.split(',')))

def get_stats(df, output_prefix):
    # reset index
    df = df.reset_index()
    ### Total reads
    df_total = df['read_id'].nunique()
    df_total_avg_length = df['read_length'].mean()
    df_total_num = len(df)
    # Unmapped reads
    unmapped = df[df['mapping_location'].str.contains('\*')]
    avg_unmapped_length = unmapped['read_length'].mean()
    unmapped_num = len(unmapped)

    # Mapped reads
    mapped = df[~df['mapping_location'].str.contains('\*')]
    avg_mapped_length = mapped['read_length'].mean()
    total_mappings = len(mapped)

    # Multi-mapped reads
    multi = mapped[mapped['mapping_location'].str.contains(',')]
    avg_multi_length = multi['read_length'].mean()
    num_multi = len(multi)

    # ASE reads
    ASE = df[~df['mapping_location'].str.contains(',') & ~df['mapping_location'].str.contains('\*')]
    av_ASE_length = ASE['read_length'].mean()
    num_ASE = len(ASE)

    #### ASE reads with mapQ > 20
    ASE_mapq = ASE[(ASE['mapq'] > 20)]
    av_ASE_mapq_length = ASE_mapq['read_length'].mean()
    num_ASE_mapq = len(ASE_mapq)
    
    ### Get stats for each haplotype
    # get uniuqe haplotypes
    haplotypes = ASE['haplotype'].unique()
    # filter ASE for each haplotype
    # Initialize an empty DataFrame
    df_all = pd.DataFrame()

    for haplotype in haplotypes:
        ASE_hap = ASE[(ASE['haplotype'] == haplotype)]
        print(ASE_hap)
        ASE_hap_avg_length = ASE_hap['read_length'].mean()
        ASE_hap_num = len(ASE_hap)
        # Create a DataFrame
        df_hap = pd.DataFrame({
            f'{output_prefix}_num': [ASE_hap_num],
            f'{output_prefix}_length': [ASE_hap_avg_length]}, index=[f'ASE_reads_{haplotype}'])
        print(df_hap)
        # Append the DataFrame to the overall DataFrame
        df_all = pd.concat([df_all, df_hap], axis=0)

    # Write the overall DataFrame to TSV
    df_all.to_csv(output_prefix + '_haplotype_summary.tsv', sep='\t')

    # Write the stats to a file 
    # Create a DataFrame
    df = pd.DataFrame({
        f'{output_prefix}_num': [df_total_num, unmapped_num, total_mappings, num_multi, num_ASE, num_ASE_mapq],
        f'{output_prefix}_length': [df_total_avg_length, avg_unmapped_length, avg_mapped_length, avg_multi_length, av_ASE_length, av_ASE_mapq_length]}, index=['Total_reads', 'Unmapped_reads', 'Mapped_reads', 'Multi_mapped_reads', 'ASE_reads', 'ASE_reads_min_mapq20'])
    print(df)
    # Write DataFrame to TSV
    df.to_csv(output_prefix + '_summary.tsv', sep='\t')

    #return total_mappings, reads_mapped, num_multi, num_ASE



def main(paf, output_prefix):
    df = paf_to_dataframe(paf) 
    df = add_haplotype_info(df, output_prefix)
    df_max = filter_max_score_keep_dups(df)
    print(df_max)

    # Get the mapping stats
    stats = get_stats(df_max, output_prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process two PAF files.')
    parser.add_argument('--paf', type=str, help='Path to the first seperate PAF file')
    parser.add_argument('--output_prefix', type=str, help='Output prefix for the summary and stats files')
    #parser.add_argument('--stats', type=str, help='Path to the stats file')
    args = parser.parse_args()
    main(args.paf, args.output_prefix)
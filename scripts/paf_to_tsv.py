import argparse
import pandas as pd

pd.options.mode.chained_assignment = None

def parse_line(line):
    fields = line.strip().split('\t')
    read_id = fields[0]
    read_length = int(fields[1])
    ref_length = int(fields[6])

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
    return read_id, read_length, mapping_location, mapq, AS_score, NM_score, ref_length

def paf_to_dataframe(file):
    data = [parse_line(line) for line in open(file, 'r')]
    df = pd.DataFrame(data, columns=['read_id','read_length', 'mapping_location', 'mapq', 'AS_score', 'NM_score', 'ref_length'])
    # return a test df
    #df = df.head(1000)
    return df

def main():
    parser = argparse.ArgumentParser(description='Convert paf file to tsv')
    parser.add_argument('paf_file', type=str, help='Input paf file')
    parser.add_argument('output_file', type=str, help='Output tsv file')
    args = parser.parse_args()

    df = paf_to_dataframe(args.paf_file)
    df.to_csv(args.output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()
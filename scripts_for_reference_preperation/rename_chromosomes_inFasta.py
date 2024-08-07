# Python
import csv
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description='Replace fasta IDs based on a mapping file.')
parser.add_argument('mapping_file', help='The mapping file.')
parser.add_argument('input_fasta', help='The input fasta file.')
parser.add_argument('output_fasta', help='The output fasta file.')
args = parser.parse_args()

# Load the mapping from the file
mapping = {}
with open(args.mapping_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        print(row)
        mapping[row[0]] = row[1]

# Read the fasta file and replace the IDs
with open(args.input_fasta, 'r') as f_in, open(args.output_fasta, 'w') as f_out:
    for line in f_in:
        if line.startswith('>'):
            id = line[0:].strip()
            if id in mapping:
                f_out.write('>' + mapping[id] + '\n')
            else:
                f_out.write(line)
        else:
            f_out.write(line)
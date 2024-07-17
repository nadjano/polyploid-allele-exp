import argparse

def add_prefix_to_fasta(input_file, output_file, prefix):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write('>' + prefix + line[1:])
            else:
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add prefix to fasta sequence IDs.')
    parser.add_argument('-i', '--input', required=True, help='Input fasta file')
    parser.add_argument('-o', '--output', required=True, help='Output fasta file')
    parser.add_argument('-p', '--prefix', required=True, help='Prefix to add to sequence IDs')
    args = parser.parse_args()

    add_prefix_to_fasta(args.input, args.output, args.prefix)
import sys
import argparse

def main():
    # Argument parsing setup
    parser = argparse.ArgumentParser(description="Process alignment and build compatibility matrix.")
    parser.add_argument('--scores', required=True, help="Path to the scores file")
    parser.add_argument('--sample', required=True, help="Sample name")
    parser.add_argument('--condition', required=True, help="conditon str")
    parser.add_argument('--output', required=True, help="Enable filtering")
    # Parse arguments
    args = parser.parse_args()

    input_file = args.scores
    output_file = args.output

    counts = {}

    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            read_id = fields[0]
            mapping_locations = []

            # Process the rest of the fields (mapping locations)
            for col in fields[1:]:
                col = col.strip()
                # Remove surrounding parentheses
                if col.startswith('(') and col.endswith(')'):
                    col = col[1:-1]
                else:
                    continue  # Skip if format is incorrect
                # Split mapping location and score
                parts = col.rsplit(',', 1)
                if len(parts) != 2:
                    continue  # Skip if format is incorrect
                mapping_location = parts[0].strip()
                try:
                    score = float(parts[1].strip())
                except ValueError:
                    continue  # Skip if score is not a number
                if score == 1.0:
                    mapping_locations.append(mapping_location)

            # Process the mapping locations with score 1.0
            if mapping_locations:
                # Sort mapping locations to ensure consistent combination keys
                mapping_locations.sort()
                # Combine mapping locations into a single key
                mapping_key = '&'.join(mapping_locations)
                counts[mapping_key] = counts.get(mapping_key, 0) + 1

    # Save the counts to the output file
    with open(output_file, 'w') as out_file:
        out_file.write(f"tname\t{args.sample}_{args.condition}\n")
        for mapping_location, count in counts.items():
            out_file.write(f"{mapping_location}\t{count}\n")

if __name__ == "__main__":
    main()
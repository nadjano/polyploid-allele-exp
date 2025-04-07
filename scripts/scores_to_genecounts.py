import sys
import argparse

def merge_dicts_to_lists(dict1, dict2):
    """
    Merge two dictionaries with overlapping keys.
    For each key in either dictionary, create a list containing values from both dictionaries.
    If a key exists in only one dictionary, use 0 for the missing value.
    """
    result = {}
    
    # Get all unique keys from both dictionaries
    all_keys = set(dict1.keys()) | set(dict2.keys())
    
    # For each key, create a list with values from both dictionaries
    for key in all_keys:
        value1 = dict1.get(key, 0)  # Use 0 if key not in dict1
        value2 = dict2.get(key, 0)  # Use 0 if key not in dict2
        result[key] = [value1, value2]
    
    return result



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

    unique_counts = {}
    multimapping_counts = {}
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

                # save the counts of the mapping locations

                # Combine mapping locations into a single key
                mapping_key = '&'.join(mapping_locations)
                counts[mapping_key] = counts.get(mapping_key, 0) + 1
                
                # if the mapping location length is > 1 the read is mulitmapping
                # split the mapping locations in unique and multimapping
                if len(mapping_locations) == 1:
                    #print(mapping_locations)
                    # Combine mapping locations into a single key
  
                    mapping_key = mapping_locations[0]
                    unique_counts[mapping_key] = unique_counts.get(mapping_key, 0) + 1
                
                elif len(mapping_locations) > 1:
                    for mapping_location in mapping_locations:
                        # Combine mapping locations into a single key
                        mapping_key = mapping_location
                        multimapping_counts[mapping_key] = multimapping_counts.get(mapping_key, 0) + 1
    
    # merge the unique and mutimapping counts to one dict
    merged_counts = merge_dicts_to_lists(unique_counts, multimapping_counts)
    # Save the counts to the output file
    with open(output_file, 'w') as out_file:
        out_file.write(f"transcript_id\tUniqueCount\tAmbigCount\n")
        for mapping_location, count in merged_counts.items():
            out_file.write(f"{mapping_location}\t{count[0]}\t{count[1]}\n")
    
    # also make stats file for number of different mapping locations, total number of reads and mulitmapping reads
    stats_file = output_file.replace('.tsv', '_stats.txt')
    with open(stats_file, 'w') as out_file:
        total_mapped_reads = sum(counts.values())  # 500
        unique_mapping_locations = len(counts)     # 5
        multimapping_reads = (sum(count for key, count in counts.items() if "&" in key)/total_mapped_reads) * 100  # 140
        unique_mapppeing_reads = (sum(count for key, count in counts.items() if "&" not in key)/total_mapped_reads) * 100
        out_file.write(f"Sample\t{args.sample}_{args.condition}\n")
        out_file.write(f"mapped_reads\t{total_mapped_reads}\n")
        out_file.write(f"unique_mapping_locations\t{unique_mapping_locations}\n")
        out_file.write(f"multimapping_reads %\t{multimapping_reads}\n")
        out_file.write(f"unique_mapping_reads %\t{unique_mapppeing_reads}\n")

if __name__ == "__main__":
    main()
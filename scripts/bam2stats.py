#!/usr/bin/env python

import sys
import pysam
import matplotlib.pyplot as plt
import argparse


def main():
    # Argument parsing setup
    parser = argparse.ArgumentParser(description="Get simple stats from bam file")
    parser.add_argument('--bam', required=True, help="Path to the bam file")
    parser.add_argument('--output', required=True, help="Sample name")
    # Parse arguments
    args = parser.parse_args()

    bam_file = args.bam
    output_prefix = args.output

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except IOError:
        print("Could not open BAM file:", bam_file)
        sys.exit(1)

    # Sets to store unique read IDs
    unique_read_ids = set()
    mapped_read_ids = set()
    unmapped_read_ids = set()

    # Dictionaries to store read lengths
    mapped_read_lengths = {}
    unmapped_read_lengths = {}

    # Lists to store read lengths for histograms
    mapped_lengths_list = []
    unmapped_lengths_list = []

    for read in bam.fetch(until_eof=True):
        read_id = read.query_name

        # Add to the set of unique read IDs
        unique_read_ids.add(read_id)
        # secondary alingments dont have a query length
        read_length = read.query_length if read.query_length != 0 else read.infer_query_length()
        
        if read.is_unmapped:
            unmapped_read_ids.add(read_id)
            # Store the read length
            unmapped_read_lengths[read_id] = read_length
            unmapped_lengths_list.append(read_length)
        else:
            mapped_read_ids.add(read_id)
            mapped_read_lengths[read_id] = read_length
            mapped_lengths_list.append(read_length)

    total_unique_reads = len(unique_read_ids)
    total_mapped_reads = len(mapped_read_ids)
    total_unmapped_reads = len(unmapped_read_ids)

    avg_mapped_length = (sum(mapped_read_lengths.values()) / len(mapped_read_lengths)
                         if mapped_read_lengths else 0)
    avg_unmapped_length = (sum(unmapped_read_lengths.values()) / len(unmapped_read_lengths)
                           if unmapped_read_lengths else 0)
    
    mapped_percentage = (total_mapped_reads / total_unique_reads) * 100
    unmapped_percentage = (total_unmapped_reads / total_unique_reads) * 100


    with open(f"{output_prefix}_stats.txt", "w") as stats_file:
        stats_file.write(f"Sample\t{output_prefix}\n")
        stats_file.write(f"Total unique read IDs\t{total_unique_reads}\n")
        stats_file.write(f"Unique mapped read IDs %\t{mapped_percentage:.2f}%\n")
        stats_file.write(f"Unique unmapped read IDs %\t{unmapped_percentage:.2f}%\n")
        stats_file.write(f"Average mapped read length\t{avg_mapped_length:.2f}\n")
        stats_file.write(f"Average unmapped read length\t{avg_unmapped_length:.2f}\n")

    # Create histograms of read lengths
    if mapped_lengths_list:
        plt.figure(figsize=(10, 6))
        plt.hist(mapped_lengths_list, bins=100, color='blue', alpha=0.7)
        plt.title('Histogram of Mapped Read Lengths')
        plt.xlabel('Read Length')
        plt.ylabel('Frequency')
        plt.xlim(0, 5000)
        plt.grid(True)
        plt.savefig(f'{output_prefix}_mapped_read_lengths_histogram.png')
        plt.close()

    # if unmapped_lengths_list:
    #     plt.figure(figsize=(10, 6))
    #     plt.hist(unmapped_lengths_list, bins=50, color='red', alpha=0.7)
    #     plt.title('Histogram of Unmapped Read Lengths')
    #     plt.xlabel('Read Length')
    #     plt.ylabel('Frequency')
    #     plt.xlim(0, max(mapped_lengths_list))
    #     plt.grid(True)
    #     plt.savefig('unmapped_read_lengths_histogram.png')
    #     plt.close()
    #     print("Unmapped read lengths histogram saved as 'unmapped_read_lengths_histogram.png'.")

if __name__ == "__main__":
    main()

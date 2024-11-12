import pysam
from tqdm import tqdm
from typing import List

class InMemoryAlignmentStore:
    def __init__(self):
        self.unique_alignments = 0
        self.groups = []

    def add_group(self, records_for_read):
        """Simulate adding a group of alignments."""
        self.groups.append(records_for_read)  # Store the group of records

    def inc_unique_alignments(self):
        """Simulate counting unique alignments."""
        self.unique_alignments += 1

class TranscriptInfo:
    """Placeholder for transcript information. Add fields as needed."""
    pass

def parse_alignments(store: InMemoryAlignmentStore, bam_file_path: str, quiet: bool = False) -> None:
    num_unmapped = 0
    prev_read = ""
    records_for_read = []

    # Open the BAM file and parse alignments
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        pb = tqdm(total=bam_file.count(until_eof=True), desc="Number of reads mapped", disable=quiet)

        for record in bam_file.fetch(until_eof=True):
            pb.update(1)
            if record.is_unmapped:
                num_unmapped += 1
                continue

            if record.query_name == prev_read:
                records_for_read.append(record)
            else:
                if prev_read:
                    store.add_group(records_for_read)
                    if len(records_for_read) == 1:
                        store.inc_unique_alignments()
                    records_for_read.clear()
                prev_read = record.query_name
                records_for_read.append(record)

        if records_for_read:
            store.add_group( records_for_read)
            if len(records_for_read) == 1:
                store.inc_unique_alignments()
        
        pb.close()

    print(f"The alignment file contained {num_unmapped} unmapped read records.")
    print(f"Unique alignments counted: {store.unique_alignments}")
    print(f"Number of groups added to store: {len(store.groups)}")

# Usage
if __name__ == "__main__":
    bam_file_path = "results/01_MINIMAP2/[SRR14993892].bam"  # Replace with the path to your BAM file
    store = InMemoryAlignmentStore()
     # Add any specific transcript information if needed

    # Parse the alignments
    parse_alignments(store, bam_file_path, quiet=False)

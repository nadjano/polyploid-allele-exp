#!/bin/bash
PROJ=/blue/mcintyre/share/potato_ASE/nf-ASE-mapping-comparison
# Define the pattern to match files
pattern="$PROJ/out_ms/scatter_plots/*/*summary_stats.tsv"

mkdir -p "$PROJ/out_ms/summary"
# Output file
output_file="$PROJ/out_ms/summary/concatenated_summary_stats.tsv"

# Initialize a flag to track the first file
first_file=true

# Loop over files matching the pattern
for file in $pattern; do
    if $first_file; then
        # Include the header for the first file
        cat "$file" > "$output_file"
        first_file=false
    else
        # Exclude the header for subsequent files
        tail -n +2 "$file" >> "$output_file"
    fi
done

echo "Concatenation complete. Output saved to $output_file."
import pandas as pd
import sys

# Get file paths from command-line arguments
file_paths = sys.argv[1:]

# Initialize an empty list to store dataframes
data_frames = []

# Loop through each file in the list of arguments and read it into a DataFrame
for file_path in file_paths:
    # Extract the sample name from the filename
    sample_name = file_path.split("/")[-1].split(".")[0]
    
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep="\t")
    
    # Rename the count column to include the sample name for clarity
    df.columns = ['tname', sample_name]
    
    # Append the DataFrame to the list
    data_frames.append(df)

# Merge all DataFrames on 'tname' column
merged_df = data_frames[0]
for df in data_frames[1:]:
    merged_df = pd.merge(merged_df, df, on='tname', how='outer')

# Save or display the merged result as needed
merged_df.to_csv("merged_counts.tsv", sep="\t", index=False)
print("Merged data saved to 'merged_counts.tsv'")
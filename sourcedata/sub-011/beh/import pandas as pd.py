import pandas as pd

# Define the input CSV file and output TSV file paths
input_csv = "sub-011_task-main_run-04_events.csv"  # Replace with your CSV file path
output_tsv = "sub-011_task-main_run-04_events.tsv"  # Replace with your desired TSV file path

# Load the CSV file
df = pd.read_csv(input_csv)

# Save it as a TSV file
df.to_csv(output_tsv, sep='\t', index=False)

print(f"Converted {input_csv} to {output_tsv}.")



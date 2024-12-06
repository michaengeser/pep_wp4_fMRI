import os
import pandas as pd

# Define the target directory
target_directory = "c:/Users/JLU-SU/OneDrive - Justus-Liebig-Universität Gießen/Dokumente/GitHub/pep_wp4_fMRI/sourcedata/sub-006/beh"

# Loop through all files in the directory
for filename in os.listdir(target_directory):
    if filename.endswith(".tsv"):  # Process only .tsv files
        file_path = os.path.join(target_directory, filename)
        
        # Open the .tsv file
        data = pd.read_csv(file_path, sep='\t')
        
        # Check and modify the 'subject' column
        if 'run' in data.columns:
            data['run'] = data['run'].apply(lambda x: '04' if x != '04' else x)
        
        # Save the modified file (overwrite)
        data.to_csv(file_path, sep='\t', index=False)
        print(f"Processed and saved {filename}.")
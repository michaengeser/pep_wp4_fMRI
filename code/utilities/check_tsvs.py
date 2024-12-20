import os
import pandas as pd

# Define the root directory
root_directory = "c:/Users/JLU-SU/OneDrive - Justus-Liebig-Universität Gießen/Dokumente/GitHub/pep_wp4_fMRI/sourcedata"

# Loop through all sub-xxx folders
for sub_folder in os.listdir(root_directory):
    if sub_folder.startswith("sub-") and os.path.isdir(os.path.join(root_directory, sub_folder)):
        subject_id = sub_folder.split("-")[1]  # Extract subject number (e.g., "001")

        # Define the 'beh' folder path
        beh_folder = os.path.join(root_directory, sub_folder, "beh")

        # Skip if 'beh' folder does not exist
        if not os.path.exists(beh_folder):
            print(f"Skipping {sub_folder}: 'beh' folder not found.")
            continue

        # Loop through .tsv files in the 'beh' folder
        for filename in os.listdir(beh_folder):
            if filename.endswith(".tsv") and "run-" in filename:
                file_path = os.path.join(beh_folder, filename)

                # Extract run number from filename
                run_number = filename.split("run-")[1][:2]  # Get the two-digit run number (e.g., "02")

                # Open the .tsv file
                data = pd.read_csv(file_path, sep='\t')

                # Check if 'subject' and 'run' columns exist
                if 'subject' in data.columns and 'run' in data.columns:
                    # Validate 'subject' column
                    if not all(data['subject'] == int(subject_id)):
                        print(f"Mismatch in {file_path}: 'subject' column does not match {subject_id}.")

                    # Validate 'run' column
                    if not all(data['run'] == int(run_number)):
                        print(f"Mismatch in {file_path}: 'run' column does not match {run_number}.")
                else:
                    print(f"Missing 'subject' or 'run' column in {file_path}.")

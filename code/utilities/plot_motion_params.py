import os
import pandas as pd
import matplotlib.pyplot as plt

# define subject, path, and runs
sub = int(input('Enter the subject number: '))
sub_id = f'sub-{sub:03d}'
sub_idx = f'sub-{sub:03d}xxxx'
path = os.path.join(
    os.path.abspath(os.path.join(os.getcwd(), '..', '..')),  
    'derivatives', sub_id, 'func'
)
runs = 13  # 1 localizer, 12 experimental runs

# load each file and plot motion parameters (translation and rotation in X/Y/Z)
for run in range(runs):
    if run == 0:  # localizer
        filename = f'rp_{sub_idx}_task-localizer_bold_00001.txt'
    else:  # experimental runs
        filename = f'rp_{sub_idx}_task-scenes_run-{run}_bold_00001.txt'

    # Full file path
    file_path = os.path.join(path, filename)

    # Load the file if it exists
    if os.path.exists(file_path):
        # open the file as space-separated values
        motion_params = pd.read_csv(file_path, sep=r'\s+', header=None)

        # plot translation (columns 1:3) and rotation (columns 4:6)
        col_ranges = [(0, 3), (3, 6)]
        titles = ['Translation', 'Rotation']

        for subplot, col_range, title in zip(range(2), col_ranges, titles):
            plt.subplot(2, 1, subplot+1)
            plt.plot(motion_params.iloc[:, col_range[0]:col_range[1]])
            plt.title(title)
            plt.legend(['X', 'Y', 'Z'])
            plt.subplots_adjust(hspace=0.5)  # add more vertical space between subplots

        # add title
        if run == 0:  # localizer
            plt.suptitle(f'Subject {sub}, Localizer Run: Translation and Rotation')
        elif run == runs - 1:  # resting state
            plt.suptitle(f'Subject {sub}, Resting State Run: Translation and Rotation')
        else:  # experimental runs
            plt.suptitle(f'Subject {sub}, Run {run}: Translation and Rotation')

        plt.show()
    else:
        print(f"File not found: {file_path}")

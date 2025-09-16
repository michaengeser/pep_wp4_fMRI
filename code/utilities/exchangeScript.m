% --- CONFIGURATION ---
backupRoot = 'D:\pep_wp4_fMRI';    % path to project backup on hard drive
laptopRoot = 'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI'; % path to project folder on laptop
filePattern = 'swrsub*bold.nii'; % files of interest

% --- FIND FILES ---
% Get list of all files starting with swsub recursively
files = dir(fullfile(backupRoot, '**', filePattern));

fprintf('Found %d files to copy.\n', numel(files));

for i = 1:numel(files)
    % Full path to the source file on the hard drive
    srcFile = fullfile(files(i).folder, files(i).name);

    % Determine relative path of the subfolder (relative to backupRoot)
    relPath = erase(files(i).folder, backupRoot);

    % Construct the matching destination folder on the laptop
    destFolder = fullfile(laptopRoot, relPath);

    % Make sure destination folder exists
    if ~exist(destFolder, 'dir')
        mkdir(destFolder);
    end

    % Construct full destination file path
    destFile = fullfile(destFolder, files(i).name);

    % Copy the file
    copyfile(srcFile, destFile);

    fprintf('Copied: %s -> %s\n', srcFile, destFile);
end

fprintf('All files copied successfully!\n');

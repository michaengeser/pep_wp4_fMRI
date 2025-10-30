% --- CONFIGURATION ---
sourceRoot = 'D:\pep_wp4_fMRI';    % i.e., path to project backup on hard drive
targetRoot = 'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\GitHub\pep_wp4_fMRI'; % i.e., path to project folder on laptop
filePattern = 'wrsub-1*'; % files of interest

% --- FIND FILES ---
% Get list of all files starting with swsub recursively
files = dir(fullfile(sourceRoot, '**', filePattern));

fprintf('Found %d files to copy.\n', numel(files));

for i = 1:numel(files)
    % Full path to the source file on the hard drive
    srcFile = fullfile(files(i).folder, files(i).name);

    % Determine relative path of the subfolder (relative to backupRoot)
    relPath = erase(files(i).folder, sourceRoot);

    % Construct the matching destination folder on the laptop
    destFolder = fullfile(targetRoot, relPath);

    % Make sure destination folder exists
    if ~exist(destFolder, 'dir')
        mkdir(destFolder);
    end

    % Construct full destination file path
    destFile = fullfile(destFolder, files(i).name);

    % Skip existing files
    if exist(destFile, 'file')
        fprintf('Skip: File %s exists already', destFile);
        disp(newline)
        continue
    end

    if contains(destFile,  {'sub-118', 'sub-125', 'sub-131'})
        fprintf('Skip %s. Excluded subject', destFile);
        disp(newline)
        continue
    end

    % Copy the file
    copyfile(srcFile, destFile);

    fprintf('Copied: %s -> %s\n', srcFile, destFile);
end

fprintf('All files copied successfully!\n');

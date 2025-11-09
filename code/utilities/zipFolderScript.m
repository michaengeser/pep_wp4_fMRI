% Define the parent folder
parentFolder = 'C:\Users\JLU-SU\OneDrive - Justus-Liebig-Universität Gießen\Dokumente\derivatives';

% Get list of all items in the parent folder
items = dir(parentFolder);

% Keep only subfolders (ignore . and ..)
subFolders = items([items.isdir] & ~ismember({items.name}, {'.', '..'}));

% Loop through each subfolder
for i = 1:numel(subFolders)
    subFolderName = subFolders(i).name;
    subFolderPath = fullfile(parentFolder, subFolderName);
    
    % Define output zip file name (in the same parent folder)
    zipFileName = fullfile(parentFolder, [subFolderName, '.zip']);
    
    fprintf('Zipping folder: %s -> %s\n', subFolderName, zipFileName);
    
    % Create the zip archive
    zip(zipFileName, subFolderPath);
end

fprintf('✅ All subfolders zipped successfully.\n');
% Due to a coding mistake the logfile of participant 001 was incorrect in
% how the trigger time stamp was stored. This was fixed with this code 
% (assisted by ChatGPT) which reconstructs the trigger time stamp using 
% the difference between the trigger time stamps and start of the first
% trial from other subjets (this may be not accurate by a few ms)


% Define the folder containing the .tsv files
% Get a list of all .tsv files in the folder
tsvFiles = dir(fullfile(pwd, '*.tsv'));

% Loop through each .tsv file
for i = 1:length(tsvFiles)
    % Construct the full path of the .tsv file
    tsvFilePath = fullfile(tsvFiles(i).folder, tsvFiles(i).name);
    
    % Read the .tsv file into a table
    T = readtable(tsvFilePath, 'FileType', 'text', 'Delimiter', '\t');
    
    % Check if 'date' and 'triggerTimeStamp' columns exist
    if ismember('date', T.Properties.VariableNames) && ismember('triggerTimeStamp', T.Properties.VariableNames)
        
        % Replace 'date' column with 'triggerDate' and copy values from 'triggerTimeStamp'
        T.Properties.VariableNames{strcmp(T.Properties.VariableNames, 'date')} = 'triggerDate';
        T.triggerDate = T.triggerTimeStamp;  % Create 'triggerDate' column
        
        % Convert 'triggerTimeStamp' from datetime to GetSecs format (using
        % other participants log files as reference)
        T.triggerTimeStamp = repmat(T.trialOnset(1) - 6.0201,...
            height(T), 1);

        % Save the modified table back to the same .tsv file
        writetable(T, tsvFilePath, 'FileType', 'text', 'Delimiter', '\t');
        
        fprintf('Processed file: %s\n', tsvFilePath);
    else
        fprintf('Skipping file (missing required columns): %s\n', tsvFilePath);
    end
end


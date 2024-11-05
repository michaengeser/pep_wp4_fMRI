%% Get Data

% Define the folder path to search for .tsv files
folderPath = pwd;  % Replace with your actual folder path

% Recursively search for all .tsv files in the folder and subfolders
tsvFiles = dir(fullfile(folderPath, '**', '*.tsv'));  % '**' searches in subfolders

% Initialize a cell array to store the data from all .tsv files
allTables = {};

% Loop through each .tsv file
for i = 1:length(tsvFiles)
    % Get the full path to the .tsv file
    tsvFilePath = fullfile(tsvFiles(i).folder, tsvFiles(i).name);
    
    % Read the .tsv file into a table
    t = readtable(tsvFilePath, 'FileType', 'text', 'Delimiter', '\t');

    % Store the table in the cell array
    if ~isempty(t)


        % make respinseKey a cell
        if ~iscell(t.responseKey)
            t.responseKey = num2cell(t.responseKey);
        end

        % add duration of run 
        runDur = t.trialEnd(end) - t.triggerTimeStamp(end) + 10;

        t.runDur = repmat(seconds(runDur), height(t), 1);
        allTables{i} = t; %#ok<SAGROW>

    end
end

% Concatenate all tables into one big table
bigTable = vertcat(allTables{:});  % Vertically concatenate all tables

%% Get Timing 
bigTable.real_stim_dur = bigTable.itiOnset - bigTable.trialOnset;
bigTable.real_stim_diff = bigTable.real_stim_dur - 0.25;
bigTable.real_trial_dur = bigTable.trialEnd - bigTable.trialOnset;
bigTable.real_trial_dur_diff = bigTable.real_trial_dur - 2.25;


% Get Timing data for each subject
timing = table('Size', [0, 8], 'VariableTypes', ...
    {'cell', 'double', 'double', 'double', 'double', 'double', 'double', 'double'},...
    'VariableNames', {'sub_num', 'mean_stim_diff', 'max_stim_diff', 'min_stim_diff',...
    'mean_trial_dur', 'max_trial_dur_diff', 'min_trial_dur_diff', 'run_dur'});

subNums = unique(bigTable.subject)';
for subNum = subNums
    timing = [timing; {num2str(subNum)}, ...
        mean(bigTable.real_stim_diff(bigTable.subject == subNum)), ...
        max(bigTable.real_stim_diff(bigTable.subject == subNum)), ...
        min(bigTable.real_stim_diff(bigTable.subject == subNum)), ...
        mean(bigTable.real_trial_dur(bigTable.subject == subNum)),...
        max(bigTable.real_trial_dur_diff(bigTable.subject == subNum)),...
        min(bigTable.real_trial_dur_diff(bigTable.subject == subNum)),...
        minutes(mean(bigTable.runDur(bigTable.subject == subNum)))];
end

% add mean and SD
timing = [timing; {'mean', ...
    mean(timing.mean_stim_diff), mean(timing.max_stim_diff),...
    mean(timing.min_stim_diff), mean(timing.mean_trial_dur),...
    mean(timing.max_trial_dur_diff), mean(timing.min_trial_dur_diff),...
    mean(timing.run_dur)}];

timing = [timing; {'sd', ...
    std(timing.mean_stim_diff(length(subNums))), std(timing.max_stim_diff(length(subNums))),...
    std(timing.min_stim_diff(length(subNums))), std(timing.mean_trial_dur(length(subNums))),...
    std(timing.max_trial_dur_diff(length(subNums))), std(timing.min_trial_dur_diff(length(subNums))),...
    std(timing.run_dur(length(subNums)))}];

% plot timing

% stimulus timing
hist_trialOnset = figure;
hist_trialOnset = histogram(bigTable.real_stim_diff, 'DisplayStyle', 'bar', 'EdgeAlpha', 0.6, 'FaceAlpha', 0.6);
title('Stimulus timing');

disp('Timin is calculated')

%% Get Accuracy 

% get table subsets
targetTable = bigTable(strcmp(bigTable.category, 'livingroom'),:);
nonTargetTable = bigTable(~strcmp(bigTable.category, 'livingroom'),:); 

% Get Timing data for each subject
accuracy = [];
subNums = unique(bigTable.subject)';
for subNum = subNums

    % Get responses for that participant
    targetResponses = targetTable.responseKey(targetTable.subject == subNum);
    targetResponses = ~strcmp(targetResponses, 'none');
    nonTargetResponses = nonTargetTable.responseKey(nonTargetTable.subject == subNum);
    nonTargetResponses = ~strcmp(nonTargetResponses, 'none');

   newTable = table;
   newTable.subNum = num2cell(subNum);
   newTable.hits =  sum(targetResponses);
   newTable.misses = sum(~targetResponses);
   newTable.hitRate = sum(targetResponses)/length(targetResponses);
   newTable.correctRejects = sum(~nonTargetResponses);
   newTable.falseAlarms = sum(nonTargetResponses);
   newTable.faRate = sum(nonTargetResponses)/length(nonTargetResponses);

   if isempty(accuracy)
       accuracy = newTable;
   else
       accuracy = [accuracy; newTable];
   end
end

% add mean and SD
newRow = cell2table({'mean', ...
    mean(accuracy.hits), mean(accuracy.misses),...
    mean(accuracy.hitRate), mean(accuracy.correctRejects),...
    mean(accuracy.falseAlarms), mean(accuracy.faRate)},...
    'VariableNames', accuracy.Properties.VariableNames);
accuracy = [accuracy; newRow];

newRow = cell2table({'sd', ...
    std(accuracy.hits(1:length(subNums))), std(accuracy.misses(1:length(subNums))),...
    std(accuracy.hitRate(1:length(subNums))), std(accuracy.correctRejects(1:length(subNums))),...
    std(accuracy.falseAlarms(1:length(subNums))), std(accuracy.faRate(1:length(subNums)))},...
     'VariableNames', accuracy.Properties.VariableNames);
accuracy = [accuracy; newRow];




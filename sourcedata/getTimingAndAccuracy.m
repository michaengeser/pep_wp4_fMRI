%% Get Data

% Define the folder path to search for .tsv files
folderPath = pwd;  % Replace with your actual folder path

% Recursively search for all .tsv files in the folder and subfolders
tsvFiles = dir(fullfile(folderPath, '**', '*.tsv'));  % '**' searches in subfolders

% Initialize a cell array to store the data from all .tsv files
allTables = {};

cfg.subNums = [101:117, 119:124, 126:130, 132, 133]; % 118, 125, and 131 did not finish the experiment
cfg = defaultCfg(cfg);

% Loop through each .tsv file
count = 0;
for i = 1:length(tsvFiles)


    if contains(tsvFiles(i).name, 'run-0') ||...
            contains(tsvFiles(i).name, 'sub-118') ||...
            contains(tsvFiles(i).name, 'sub-125') ||...
            contains(tsvFiles(i).name, 'sub-131')
        continue
    end
    count = count + 1;


    % Get the full path to the .tsv file
    tsvFilePath = fullfile(tsvFiles(i).folder, tsvFiles(i).name);
    
    % Read the .tsv file into a table
    t = readtable(tsvFilePath, 'FileType', 'text', 'Delimiter', '\t');

    % Store the table in the cell array
    if ~isempty(t)

        % make responseKey a cell
        if ~iscell(t.responseKey)
            t.responseKey = num2cell(t.responseKey);
        end

        % add duration of run 
        runDur = t.trialEnd(end) - t.triggerTimeStamp(end) + 10;

        t.runDur = repmat(seconds(runDur), height(t), 1);
        allTables{count} = t; %#ok<SAGROW>

    end
end

% Concatenate all tables into one big table
bigTable = vertcat(allTables{:});  % Vertically concatenate all tables
bigTable = bigTable(bigTable.subject ~= 131, :);

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

% Get Timing data for each subject
accuracy = [];
verticalTable = [];
subNums = unique(bigTable.subject)';
for subNum = subNums

    % get table subsets
    targetTable = bigTable(~isnan(bigTable.accuracy) & bigTable.subject == subNum,:);
    subjecTable = bigTable(bigTable.subject == subNum,:);

    % make binomial test against chance
    x = sum(targetTable.accuracy);           
    n = length(targetTable.accuracy);           
    pNull = 0.5;      
    pVal = 1 - binocdf(x - 1, n, pNull);  % subtract 1 because binocdf is P(X ≤ x)


    newTable = table;
    newTable.subNum = num2cell(subNum);
    newTable.numCorrect = sum(targetTable.accuracy == 1);
    newTable.numIncorrect = sum(targetTable.accuracy == 0);
    newTable.meanRT = mean(targetTable.responseTime - targetTable.trialEnd);
    newTable.meanAccuracy = mean(targetTable.accuracy);
    newTable.pVal = pVal;

    if isempty(accuracy)
        accuracy = newTable;
        verticalTable = subjecTable;
    else
        accuracy = [accuracy; newTable];
        %
        newNames = strcat(subjecTable.Properties.VariableNames, "_", string(repmat(subNum, 1, width(subjecTable))));
        subjecTable.Properties.VariableNames = newNames;
        %verticalTable = [verticalTable, subjecTable];
    end
end

% add mean and SD
newRow = cell2table({'mean', ...
    mean(accuracy.numCorrect), mean(accuracy.numIncorrect),...
    mean(accuracy.meanRT), mean(accuracy.meanAccuracy), mean(accuracy.pVal)},...
    'VariableNames', accuracy.Properties.VariableNames);
accuracy = [accuracy; newRow];

newRow = cell2table({'sd', ...
    std(accuracy.numCorrect(1:length(subNums))), std(accuracy.numIncorrect(1:length(subNums))),...
    std(accuracy.meanRT(1:length(subNums))), std(accuracy.meanAccuracy(1:length(subNums))), std(accuracy.pVal(1:length(subNums)))},...
     'VariableNames', accuracy.Properties.VariableNames);
accuracy = [accuracy; newRow];



%% other stuff

% verticalTable2 = verticalTable(~isnan(verticalTable.accuracy_102), :);
% verticalTable3 = verticalTable(: ,contains(verticalTable.Properties.VariableNames, 'acc'));
% cor_rdm = corr(table2array(verticalTable3), 'rows', 'pairwise');
% 
% % get target structure
% targetDir = fullfile('pwd', '..', 'code', 'utilities', 'targets.mat');
% load(targetDir)
% 
% warning('off')
% endTable = [];
% for i = 1:12
% 
%     newTable = table;
%     for ii = 1:15
%         newTable.imageName(ii) = targetStruct(i).imgName(ii);
%         newTable.trialNum(ii) = targetStruct(i).trialNum(ii);
%         newTable.targetName(ii) = targetStruct(i).targetName(ii);
%         newTable.targetPresent(ii) = targetStruct(i).targetPresent(ii);
%         newTable.run(ii) = i;
%     end
% 
%     if isempty(endTable)
%         endTable = newTable;
%     else
%         endTable = [endTable; newTable];
%     end
% end
% warning('on')


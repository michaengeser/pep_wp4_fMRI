load(fullfile(pwd, 'utilities', 'targets.mat'), 'targetStruct')

% Initialize storage for unique strings and their counts
allNames = [];  % Store all target names
runIndices = [];  % Store corresponding run indices
targetPresentCounts = [];  % Store number of times target was present

% Loop through all runs
for run = 1:12
    % Extract target names and targetPresent values for current run
    targetNames = targetStruct(run).targetName;  % Cell array of strings
    targetPresent = targetStruct(run).targetPresent;  % Logical array

    % Append target names and their corresponding run indices
    allNames = [allNames; targetNames'];
    runIndices = [runIndices; repmat(run, length(targetNames), 1)];

    % Count occurrences where target is present (true)
    targetPresentCounts = [targetPresentCounts; targetNames(targetPresent)'];
end

% Get unique target names
[uniqueNames, ~, idx] = unique(allNames);

% Initialize result storage
numUnique = length(uniqueNames);
occurrences = zeros(numUnique, 1);  % Count occurrences
presentCounts = zeros(numUnique, 1);  % Count true occurrences
runsList = cell(numUnique, 1);  % Store run indices

% Process each unique name
for i = 1:numUnique
    % Get indices of occurrences for this unique name
    nameMask = strcmp(allNames, uniqueNames{i});

    % Count occurrences
    occurrences(i) = sum(nameMask);

    % Get unique runs where this name appeared
    runsList{i} = unique(runIndices(nameMask))';  % Transpose for row format

    % Count occurrences where targetPresent was true
    nameMask = strcmp(targetPresentCounts, uniqueNames{i});
    presentCounts(i) = sum(nameMask);

end

% Create results table
resultsTable = table(uniqueNames, occurrences, presentCounts, runsList, ...
    'VariableNames', {'TargetName', 'TotalOccurrences', 'PresentCount', 'Runs'});

% Display table
disp(resultsTable);

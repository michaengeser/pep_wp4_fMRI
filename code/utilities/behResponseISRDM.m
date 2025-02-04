function d = behResponseISRDM(cfg, d)

% evaluate input
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85;end
if ~isfield(cfg, 'nVols'); cfg.nVols = 152;end
if ~isfield(cfg, 'nTrials'); cfg.nTrials = 100;end
if ~isfield(cfg, 'nTargets'); cfg.nTargets = 16; end

% init variables
subs = cfg.subNums;

%% Define subjects and main path
mainPath = fullfile(pwd, '..');
behPath = fullfile(mainPath, 'sourcedata');

% initialize subjects response table
nTrialsPerBlock = (cfg.nTrials + cfg.nTargets) / 2;
bathroomDataMat = nan(cfg.nRuns * nTrialsPerBlock, cfg.n);
kitchenDataMat = bathroomDataMat;


%% get data

for iSub = 1:length(subs)


    % runtime control
    disp(['Subject: ', num2str(subs(iSub))])

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    %
    subFolder = fullfile(behPath, subID, 'beh');
    files = dir(fullfile(subFolder, '*.tsv'));
    files = struct2table(files);

    %% Loop through runs

    for iRun = 1:cfg.nRuns
        
        % for subjects 1 in run 10 the kitchen block was not shown
        % correctly, skip that
        if subs(iSub) == 1 && iRun == 10
            continue
        end

        % get data
        fileData = readtable(char(fullfile(files.folder(iRun), files.name(iRun))),...
            'FileType', 'text', 'Delimiter', '\t');

        if strcmp(fileData.category{1}, 'bathroom')
            bathroomTable = fileData(1:nTrialsPerBlock, :);
            kitchenTable = fileData(nTrialsPerBlock + 1:end, :);
        else
            kitchenTable = fileData(1:nTrialsPerBlock, :);
            bathroomTable = fileData(nTrialsPerBlock + 1:end, :);
        end

        % store response vector in matrix
        idx = 1 + (iRun-1) * nTrialsPerBlock : iRun * nTrialsPerBlock;
        kitchenDataMat(idx, iSub) = ~isnan(kitchenTable.responseTime)';
        bathroomDataMat(idx, iSub) = ~isnan(bathroomTable.responseTime)';
       
    end % runs
end % subjects

% get IS RDM
[~, d.behResponseRDM.bathroom.IS_RDM.RDM, ~] = make_RDM(bathroomDataMat, cfg);
d.behResponseRDM.bathroom.IS_RDM.name = 'IS_RDM_beh_responses';
d.behResponseRDM.bathroom.IS_RDM.color = [0, 0, 0];
[~, d.behResponseRDM.kitchen.IS_RDM.RDM, ~] = make_RDM(kitchenDataMat, cfg);
d.behResponseRDM.kitchen.IS_RDM.name = 'IS_RDM_beh_responses';
d.behResponseRDM.kitchen.IS_RDM.color = [0, 0, 0];

d.behResponseRDM.bathroom.subMatrix = bathroomDataMat;
d.behResponseRDM.kitchen.subMatrix = kitchenDataMat;

end
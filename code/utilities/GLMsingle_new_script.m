

nRuns = 10;
tr = 1.85;
nTrials = 100;
nVols = 152;
stimdur = 0.25;



%% make multiple condition files
sortRows = true;
includeTargets = false;
create_mcf_func(subs, sortRows, includeTargets)


spm_jobman('initcfg');

%% Define subjects and main path
mainPath = fullfile(pwd, '..');
fmriPath = fullfile(mainPath, 'sourcedata');



%% get data

for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    % make `derivatives` sub-directory if it doesn't exist yet
    if ~exist(fullfile(mainPath, 'derivatives', subID), 'dir')
        mkdir(fullfile(mainPath, 'derivatives', subID));
    end

    % make make output fodler doesn't exist yet
    outputdir = fullfile(mainPath, 'derivatives', subID, 'GLMsingleEstimates');
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end

    %% Loop through runs

    % for each run, get each scan's `.nii` file
    data_new = cell(1, nRuns);
    design = cell(1, nRuns);
    trialIDs = [];
    %motion = cell(1, nRuns); % not  necessary

    for iRun = 1:nRuns

        % get functional files
        funcFile = fullfile(mainPath, 'derivatives', subID, 'func', ...
            sprintf('wr%s%s_task-scenes_run-%s_bold.nii',...
            subID, 'xxxx', num2str(iRun)));

        if ~exist(funcFile, 'file')
            error('make 4D file')
        end

        % load data
        data_new{iRun} = load_untouch_nii(funcFile);

        % get `mcf` file
        mcf = fullfile(fmriPath, subID, 'beh', 'onsets', ...
            sprintf('mcf_%s_run-%s.mat', subID, num2str(iRun)));
        mcf_file = load(mcf);

        % add trial IDs
        trialIDs = [trialIDs; cell2mat(mcf_file.trialIDs),...
             ones(size(mcf_file.trialIDs))*iRun];


        % check if number of conditions match onsets
        if ~numel(mcf_file.onsets) == nTrials
            warning('Onsets does not match number of conditions')
        end

        % make design matrix
        design{iRun} = zeros(nVols,nTrials);

        % set onset in correct TR
        for iTrials=1:nTrials     %
            design{iRun}(round(cell2mat(mcf_file.onsets(iTrials)) / tr) + 1,...
                iTrials) = 1;
        end

        %     % get motion regressors
        %     moRegs = fullfile(mainPath, 'derivatives', subID, 'func',...
        %         sprintf('rp_%s%s_task-scenes_run-%s_bold_00001.txt',...
        %         subID, 'xxxx', num2str(iRun)));
        %     motion{iRun} = load(moRegs);

    end

    % set options
    opt = struct('wantmemoryoutputs',[0 0 0 1]); % only type d model should be written to memory
    opt.wantfileoutputs = [0 0 0 1];
    opt.extraregressors = motion;

    % run GLM single
    [results designSINGLE] = GLMestimatesingletrial(design,data,stimdur,tr,outputdir,opt);

    %% get reliability

    % Create output variable for reliability values
    pairs = nchoosek(1:10,2);
    vox_reliabilities = cell(1,length(pairs));

    % For each pair of trials
    for p = 1 : length(pairs)

        % Get the GLM betas
        betas = models.(model_names{p}).modelmd(:,:,:,repindices);  % use indexing to pull out the trials we want
        betas_reshaped = reshape(betas,size(betas,1),size(betas,2),size(betas,3),2,[]);  % reshape to X x Y x Z x 2 x CONDITIONS

        % compute reliabilities using an efficient (vectorized) utility
        % function
        vox_reliabilities{p} = calccorrelation(betas_reshaped(:,:,:,1,:),betas_reshaped(:,:,:,2,:),5);

        % Note that calccorrelation.m is a utility function that computes
        % correlations in a vectorized fashion (for optimal speed).

    end


end


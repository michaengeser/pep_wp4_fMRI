function ds_out = loadCosmoDataset(cfg, subID, mask_fn)

if strcmp(cfg.map, 't')
    % Initialize dataset cell array
    nTrials = cfg.nTrials;
    nTargetIDperRun = nTrials/2;
    datasets = cell(nTargetIDperRun, cfg.nRuns);

    % Loop over runs to load data
    if isempty(gcp('nocreate'))
        parpool(8);
    end
    parfor iRun = 1:cfg.nRuns

        % get dir for folder
        con_dir = fullfile(pwd, '..', 'derivatives', subID, 'exp_glm1_norm', ...
            sprintf('run-%02d',iRun));

        % load trial IDs
        beh_dir = fullfile(pwd, '..', 'sourcedata', subID, 'beh', ...
            'onsets', sprintf('mcf_%s_run-%d', subID, iRun));
        targetID = load(beh_dir, 'trialIDs')
        runTrials = unique(cell2mat(targetID.trialIDs))';

        for trial = 1:nTargetIDperRun
            trialID = runTrials(trial);

            % Load the contrast image for this run and trial
            con_file = fullfile(con_dir, sprintf('con_%04d.nii',trial));

            % Convert to CoSMoMVPA dataset
            ds = cosmo_fmri_dataset(con_file, ...
                'mask', mask_fn, ... % Set brain mask
                'targets', trialID, ... % Set target ID
                'chunks', ceil(iRun/2));     % Set run identifiers

            % Store dataset
            datasets{trial, iRun} = ds;
        end
    end

    % Combine all runs into a single dataset
    ds_per_run = cosmo_stack(datasets);

elseif strcmp(cfg.map, 'b')

    % get path
    betaPath = fullfile(pwd, '..', 'derivatives', subID, 'GLMsingle_betas', ...
        'beta_sorted.nii');

    % load beta map
    ds_per_run = cosmo_fmri_dataset(betaPath, ...
        'mask', mask_fn); % Set brain mask

    % add chunks and targets
    nSamples = height(ds_per_run.samples);
    nSamplesPerRun = nSamples/cfg.nRuns;
    ds_per_run.sa.targets = repmat(1:cfg.nTrials, 1, cfg.nRuns)';
    preChunks = repmat(1:cfg.nRuns, nSamplesPerRun, 1);
    ds_per_run.sa.chunks = reshape(preChunks,[],1);


else
    error('map not defined')
end

% remove constant features
if iscell(ds_per_run.sa.targets)
    ds_per_run.sa.targets = cell2mat(ds_per_run.sa.targets);
end
ds_out=cosmo_remove_useless_data(ds_per_run);

% check the dataset
cosmo_check_dataset(ds_out);
end
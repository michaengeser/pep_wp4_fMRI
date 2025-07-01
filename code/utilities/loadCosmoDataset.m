function ds_out = loadCosmoDataset(cfg, subID, mask_label)

% evaluate input
if ~isfield(cfg, 'pca'); cfg.pca = true; end

% get mask directory
if contains(mask_label, 'LPFC')
    funcROIname = [mask_label(2:5), '_funcROI.nii'];
else
    funcROIname = [mask_label(2:4), '_funcROI.nii'];
end
funcROIdir = fullfile(pwd, '..', 'MNI_ROIs', 'func_ROIs', subID, funcROIname);

% if functional defined ROI available take that
if exist(funcROIdir, 'file')
    mask_fn = funcROIdir;
else
    mask_fn = fullfile(pwd, '..', 'MNI_ROIs', [char(mask_label)]);
end

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
    betaPath = fullfile(pwd, '..', 'derivatives', subID, 'GLMsingleEstimates', ...
        'GLMsingle_betas.nii');

    % load beta map
    ds_per_run = cosmo_fmri_dataset(betaPath, ...
        'mask', mask_fn); % Set brain mask

    % add chunks and targets
    ids = load(fullfile(pwd, '..', 'derivatives', subID, 'GLMsingleEstimates', ...
        'trialIDs.mat'));
    ds_per_run.sa.targets = ids.trialIDs(:, 1);
    ds_per_run.sa.chunks = ceil(ids.trialIDs(:, 2)/2);

else
    error('map not defined')
end

% remove constant features
if iscell(ds_per_run.sa.targets)
    ds_per_run.sa.targets = cell2mat(ds_per_run.sa.targets);
end
ds_out=cosmo_remove_useless_data(ds_per_run);

% reduce number of features using PCA
if cfg.pca
    [ds_out, pca_params] = cosmo_map_pca(ds_out, 'max_feature_count', 11000, 'pca_explained_ratio', 0.95);
    ds_out.sa = ds_per_run.sa;
end

% check the dataset
cosmo_check_dataset(ds_out);
end
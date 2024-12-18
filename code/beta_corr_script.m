% Initialize dataset cell array

mask_fn=fullfile(pwd, '..', 'MNI_ROIs', 'wlV1.nii');

nRuns = 10;

for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    % Loop over runs to load data
    counter = 0;
    datasets = {};
    for iRun = 1:nRuns

        % get dir for folder
        con_dir = fullfile(pwd, '..', 'derivatives', subID, 'exp_glm1_norm', ...
            sprintf('run-%02d',iRun));

        % get files
        files = dir(fullfile(con_dir, 'con*.nii'));

        % get behavioral log file
        runID = sprintf('run-%02d', iRun);
        logFile = readtable(fullfile(pwd, '..', 'sourcedata', subID, 'beh', ...
            sprintf('%s_task-main_%s_events.tsv', subID, runID)),...
            'FileType', 'text', 'Delimiter', '\t');

        logFile = logFile(~strcmp(logFile.category, 'livingroom'),:);

        for iFile = 1:height(logFile)
            counter = counter + 1;

            % Load the contrast image for this run and trial
            con_file = fullfile(con_dir, files(iFile).name);

            % Convert to CoSMoMVPA dataset
            ds = cosmo_fmri_dataset(con_file, ...
                'mask', mask_fn, ... % Set brain mask
                'targets', logFile.texture(iFile) - 10, ...
                'chunks', logFile.run(iFile));     % Set run identifiers

            % Store dataset
            datasets{counter} = ds;
        end
    end

    % Combine all runs into a single dataset
    ds_per_run = cosmo_stack(datasets);

    % remove constant features
    ds_per_run = cosmo_remove_useless_data(ds_per_run);


    % Combine chunks and targets into a sorting matrix
    sort_matrix = [ds_per_run.sa.chunks(:), ds_per_run.sa.targets(:)];
    [~, sort_idx] = sortrows(sort_matrix);
    ds_per_run.samples = ds_per_run.samples(sort_idx, :);
    ds_per_run.sa.chunks = ds_per_run.sa.chunks(sort_idx);
    ds_per_run.sa.targets = ds_per_run.sa.targets(sort_idx);


    rdm_corr = corr(ds_per_run.samples', 'type', 'Spearman');

    subID2 = strrep(subID, '-', '');
    RDMstruct.(subID2) = rdm_corr;

    figure; imagesc(rdm_corr, [-1,1]); colorbar
 title('Pairwise Correlation in COMSMO');
end





%----------------------------------------------------------



betaPath = fullfile(pwd, '..', 'derivatives', subID, 'GLMsingle_betas', ...
    'beta_sorted.nii');
ds_per_run = cosmo_fmri_dataset(betaPath, ...
    'mask', mask_fn); % Set brain mask


% add chunks and targets
nSamples = height(ds_per_run.samples);
nSamplesPerRun = nSamples/nRuns;
ds_per_run.sa.targets = repmat(1:nTrials, 1, nRuns)';
preChunks = repmat(1:nRuns, nSamplesPerRun, 1);
ds_per_run.sa.chunks = reshape(preChunks,[],1);

nFeartures = size(ds_per_run.samples, 2);
features_all = zeros(nTrials, nFeartures, nRuns);

% Create a tiled layout with 10 tiles (1 for each run)
figure
tiledlayout(2, 5); % Adjust rows and columns as necessary for better visualization

% Loop through each run to compute and display the RDM
for run = 1:nRuns

    % all features
    features_all(:, :, run) = ds_per_run.samples(ds_per_run.sa.chunks' == run, :);

    % Extract features for the current run
    features_run = ds_per_run.samples(ds_per_run.sa.chunks' == run, :);

    % Compute the RDM (correlation matrix) for the current run
    rdm_corr_run = corr(features_run', 'type', 'Spearman');

    % Plot the RDM in the corresponding tile
    nexttile;
    imagesc(rdm_corr_run, [-1, 1]);
    title(['Run ', num2str(run)]);
    colorbar;
    axis square; % Make the plots square
end
sgtitle('RDMs for Each Run'); % Set a global title for the tiled plot


features_mean = mean(features_all, 3);

rmd_corr = corr(features_mean', 'type', 'Spearman');
figure; imagesc(rmd_corr, [-1,1]);


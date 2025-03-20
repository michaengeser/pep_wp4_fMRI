function res = PairwiseROIsLDAclassifier(cfg)


% Define classifiers
classifier = @cosmo_classify_lda;

% get number of ROI masks
nmasks=numel(cfg.rois);

% loop through subjects
for iSub = 1:length(cfg.subNums)

    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

    % progress report
    disp(['Starting pairwise ROI decoding for subject ',  num2str(cfg.subNums(iSub)), ' on ',...
        cfg.map, '-map']);


    for j=1:nmasks

        mask_label=cfg.rois{j};
        mask_fn=fullfile(pwd, '..', 'MNI_ROIs', [char(mask_label)]);

        disp(['Using mask ',  mask_label]);
        disp(char(datetime))

        if strcmp(cfg.map, 't')
            % Initialize dataset cell array
            datasets = cell(cfg.nTrials, cfg.nRuns);

            % Loop over runs to load data
            nTrials = cfg.nTrials;
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

                for trial = 1:nTrials

                    % Load the contrast image for this run and trial
                    con_file = fullfile(con_dir, sprintf('con_%04d.nii',trial));

                    % Convert to CoSMoMVPA dataset
                    ds = cosmo_fmri_dataset(con_file, ...
                        'mask', mask_fn, ... % Set brain mask
                        'targets', targetID.trialIDs(trial), ... % Set condition labels (1 = bathroom)
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
        ds=cosmo_remove_useless_data(ds_per_run);

        % reduce number of features using PCA
        [ds, pca_params] = cosmo_map_pca(ds, 'max_feature_count', 10000, 'pca_explained_ratio', 0.99);
        ds.sa = ds_per_run.sa;
        cosmo_check_dataset(ds);

        %% Pairwise decoding
        % Initialize RDM
        rdm = zeros(cfg.nTrials, cfg.nTrials);
        disp('')

        % define options
        opt.max_feature_count = 5000;

        if isempty(gcp('nocreate'))
            parpool(8);
        end
        parfor s1 = 1:nTrials
            for s2 = 1:nTrials
                if ~(s2 > s1)
                    continue
                end

                % Subset data for the two stimuli
                ds_stim = cosmo_slice(ds, ds.sa.targets == s1 | ds.sa.targets == s2);

                % Rename target
                ds_stim.sa.targets = (ds_stim.sa.targets == s1) + 1;

                % Define partitions
                partitions=cosmo_nfold_partitioner(ds_stim);

                % get predictions for each fold
                [~ ,accuracy] = cosmo_crossvalidate(ds_stim, classifier, partitions, opt);

                % Store the accuracy
                rdm(s2, s1) = accuracy;                
            end
        end

        % make rdm symetric
        rdm = squareform(squareform(rdm));

        % Save the RDM
        mask_label_short = split(mask_label, '.');
        mask_label_short = mask_label_short{1};
        subID2 = strrep(subID, '-', '');

        res.(subID2).(mask_label_short).rdm = rdm;
        res.(subID2).(mask_label_short).mean_accuracy = mean(squareform(rdm));
    end
end

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

%% ploting

% Get list of subjects and ROIs
subjects = fieldnames(res);
masks = fieldnames(res.(subjects{1}));

% Initialize data storage
num_subjects = numel(subjects);
num_rois = numel(masks);
all_data = nan(num_subjects, num_rois);

% Collect data
for i_sub = 1:num_subjects
    subID = subjects{i_sub};
    for i_roi = 1:num_rois
        mask_label = masks{i_roi};
        all_data(i_sub, i_roi) = res.(subID).(mask_label).mean_accuracy;
    end
end

% Compute mean and standard deviation for each ROI
mean_data = mean(all_data, 1);
std_data = std(all_data, 0, 1);

% Create bar plot
figure;
hold on;

% Bar plot with error bars
bar_handle = bar(mean_data, 'FaceColor', 'flat');
errorbar(1:num_rois, mean_data, std_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add horizontal line at chance level
yline(0.5, '--r', 'Chance Level', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');

% Add jittered individual points
rng(0); % For reproducible jitter
jitter_amount = 0.1; % Adjust jitter spread
for i_roi = 1:num_rois
    x_jitter = i_roi + (rand(num_subjects, 1) - 0.5) * jitter_amount;
    scatter(x_jitter, all_data(:, i_roi), 30, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
    for i_sub = 1:num_subjects
        text(x_jitter(i_sub) + 0.04, all_data(i_sub, i_roi), subjects{i_sub}, 'FontSize', 8);
    end
end

% Customize plot
xticks(1:num_rois);
xticklabels(masks);
xlabel('ROI');
ylabel('Accuracy');
ylim([0.4,0.6])
title('Mean Pairwise Decoding Results Across ROIs');

hold off;



end
function res = ROIsLDAclassifier(subs, nRuns, nTrials, map, rois)


% Define classifiers
classifier = @cosmo_classify_lda;

% get number of ROI masks
nmasks=numel(rois);

% loop through subjects
for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    disp(['Start ROI decoding for subject ',  num2str(subs(iSub)), ' on ',...
        map, '-map']);


    for j=1:nmasks

        mask_label=rois{j};
        mask_fn=fullfile(pwd, '..', 'MNI_ROIs', [char(mask_label), '.img']);

        if strcmp(map, 't')
            % Initialize dataset cell array
            datasets = cell(numel(nRuns), 1);

            % Loop over runs to load data
            counter = 0;
            for run = 1:nRuns

                % get dir for folder
                con_dir = fullfile(pwd, '..', 'derivatives', subID, 'exp_glm1_norm', ...
                    sprintf('run-%02d',run));

                for trial = 1:nTrials
                    counter = counter + 1;

                    % Load the contrast image for this run and trial
                    con_file = fullfile(con_dir, sprintf('con_%04d.nii',trial));

                    % Convert to CoSMoMVPA dataset
                    if trial <= nTrials/2; target = 1; else; target = 2;end
                    ds = cosmo_fmri_dataset(con_file, ...
                        'mask', mask_fn, ... % Set brain mask
                        'targets', target, ... % Set condition labels (1 = bathroom)
                        'chunks', run);     % Set run identifiers

                    % Store dataset
                    datasets{counter} = ds;
                end
            end

            % Combine all runs into a single dataset
            ds_per_run = cosmo_stack(datasets);

        elseif strcmp(map, 'b')

            % get path
            betaPath = fullfile(pwd, '..', 'derivatives', subID, 'GLMsingle_betas', ...
                'beta.nii');

            % load beta map
            ds_per_run = cosmo_fmri_dataset(betaPath, ...
                'mask', mask_fn); % Set brain mask

            % add chunks and targets
            nSamples = height(ds_per_run.samples);
            nSamplesPerRun = nSamples/nRuns;
            ds_per_run.sa.targets = repmat([ones(1,nTrials/2), ones(1,nTrials/2)*2,...
                ones(1,nSamplesPerRun-nTrials)*3], 1, nRuns)';
            preChunks = repmat(1:nRuns, nSamplesPerRun, 1);
            ds_per_run.sa.chunks = reshape(preChunks,[],1);

            % remove living room trials trials (target == 3)
            ds_per_run = cosmo_slice(ds_per_run, (ds_per_run.sa.targets ~= 3), 1);

        else
            error('map not defined')
        end

        % remove constant features
        ds=cosmo_remove_useless_data(ds_per_run);

        % print dataset
        fprintf('Dataset input:\n');
        cosmo_disp(ds);

        % Define partitions
        partitions=cosmo_nfold_partitioner(ds);

        % print dataset
        fprintf('Partitions:\n');
        cosmo_disp(partitions);

        % get predictions for each fold
        [pred,accuracy] = cosmo_crossvalidate(ds, classifier, partitions);

        % get confusion matrix for each fold
        confusion_matrix_folds=cosmo_confusion_matrix(ds.sa.targets,pred);

        % sum confusion for each ground-truth target and prediction,
        % resulting in an nclasses x nclasses matrix
        confusion_matrix = sum(confusion_matrix_folds,3);

        % store results
        subID2 = strrep(subID, '-', '');
        res.(subID2).(mask_label).confusion_matrix = confusion_matrix / (nRuns*nTrials);
        res.(subID2).(mask_label).accuracy = accuracy;

    end
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
        all_data(i_sub, i_roi) = res.(subID).(mask_label).accuracy;
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
        text(x_jitter(i_sub) + 0.03, all_data(i_sub, i_roi), subjects{i_sub}, 'FontSize', 8);
    end
end

% Customize plot
xticks(1:num_rois);
xticklabels(masks);
xlabel('ROI');
ylabel('Accuracy');
ylim([0.4,0.6])
title('Decoding Results Across ROIs');

hold off;



end
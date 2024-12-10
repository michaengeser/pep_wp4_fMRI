function res = CrossValdidatedRDM(subs, nRuns, nTrials, map, rois)


%% get data 

% get number of ROI masks
nmasks=numel(rois);

% loop through subjects
for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    % progress report
    disp(['Starting cross validated RDMg for subject ',  num2str(subs(iSub)), ' on ',...
        map, '-map']);


    for j=1:nmasks

        mask_label=rois{j};
        mask_fn=fullfile(pwd, '..', 'MNI_ROIs', [char(mask_label)]);

        disp(['Using mask ',  mask_label]);
        disp(char(datetime))

        if strcmp(map, 't')
            % Initialize dataset cell array
            datasets = cell(numel(nRuns), 1);

            % Loop over runs to load data
            counter = 0;
            for run = 1:nRuns
                disp(num2str(run));
                % get dir for folder
                con_dir = fullfile(pwd, '..', 'derivatives', subID, 'exp_glm1_norm_with_target/', ...
                    sprintf('run-%02d',run));

                for trial = 1:nTrials
                    counter = counter + 1;

                    % Load the contrast image for this run and trial
                    con_file = fullfile(con_dir, sprintf('con_%04d.nii',trial));

                    % Convert to CoSMoMVPA dataset
                    ds = cosmo_fmri_dataset(con_file, ...
                        'mask', mask_fn, ... % Set brain mask
                        'targets', trial, ... % Set condition labels (1 = bathroom)
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
            ds_per_run.sa.targets = repmat(1:nTrials, 1, nRuns)';
            preChunks = repmat(1:nRuns, nSamplesPerRun, 1);
            ds_per_run.sa.chunks = reshape(preChunks,[],1);


        else
            error('map not defined')
        end

        % remove living room trials trials (target > 100)
        ds_per_run = cosmo_slice(ds_per_run, (ds_per_run.sa.targets <= 100), 1);

        % remove constant features
        ds=cosmo_remove_useless_data(ds_per_run);

        %% Crossvalidated RDM 

        % Initialize RDM
        rdm = zeros(nTrials, nTrials);

        % get combination to split 10 into 2 halfs
        splits = nchoosek(1:10, 5);
        splits = splits(1:height(splits)/2, :);

        % just do even/odd for now
        %splits = [1,3,5,7,9];

        disp('')
        for s1 = 1:nTrials
            for s2 = s1+1:nTrials

                % Subset data for the two stimuli
                ds_stim = cosmo_slice(ds, ds.sa.targets == s1 | ds.sa.targets == s2);

                % Rename target
                ds_stim.sa.targets = (ds_stim.sa.targets == s1) + 1;

                % loop over splits
                crossValdidatedRs = zeros(1, height(splits));
                for iSplit = 1:height(splits)
                    split = splits(iSplit, :);

                    % split data
                    dataSplit1Cond1 = cosmo_slice(ds_stim, ds_stim.sa.targets == 1 &...
                        ismember(ds_stim.sa.chunks, split));
                    dataSplit1Cond2 = cosmo_slice(ds_stim, ds_stim.sa.targets == 2 &...
                        ismember(ds_stim.sa.chunks, split));
                    dataSplit2Cond1 = cosmo_slice(ds_stim, ds_stim.sa.targets == 1 &...
                        ~ismember(ds_stim.sa.chunks, split));
                    dataSplit2Cond2 = cosmo_slice(ds_stim, ds_stim.sa.targets == 2 &...
                        ~ismember(ds_stim.sa.chunks, split));

                    % take mean for each split
                    meanData = zeros(width(dataSplit1Cond1.samples), 4);
                    meanData(:, 1) = mean(dataSplit1Cond1.samples)';
                    meanData(:, 2) = mean(dataSplit1Cond2.samples)';
                    meanData(:, 3) = mean(dataSplit2Cond1.samples)';
                    meanData(:, 4) = mean(dataSplit2Cond2.samples)';

                    % correlate them
                    rvals = corr(meanData, 'type', 'Spearman');

                    % get crossvalidation R
                    within = mean([rvals(1, 3), rvals(2, 4)]); % same cond, other split
                    between = mean([rvals(1, 4), rvals(2, 3)]); % other cond, other split
                    crossValdidatedRs(iSplit) =  between - within;

                end 


                % Store the correlation
                rdm(s1, s2) = mean(crossValdidatedRs);
                rdm(s2, s1) = rdm(s1, s2); % Symmetric matrix
            end
        end

        % Save the RDM
        mask_label_short = split(mask_label, '.');
        mask_label_short = mask_label_short{1};
        subID2 = strrep(subID, '-', '');

        res.(subID2).(mask_label_short).rdm = rdm;
        res.(subID2).(mask_label_short).mean_cor = mean(squareform(rdm));
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
        all_data(i_sub, i_roi) = res.(subID).(mask_label).mean_cor;
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
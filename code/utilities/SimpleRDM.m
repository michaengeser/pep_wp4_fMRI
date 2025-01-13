function res = SimpleRDM(subs, nTrials, map, rois)

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
    disp(['Starting spm RDM for subject ',  num2str(subs(iSub)), ' on ',...
        map, '-map']);

    % loop through ROIs
    for j=1:nmasks

        mask_label=rois{j};
        mask_fn=fullfile(pwd, '..', 'MNI_ROIs', [char(mask_label)]);
        mask_label_short = split(mask_label, '.');
        mask_label_short = mask_label_short{1};


        disp(['Using mask ',  mask_label]);
        disp(char(datetime))

        if strcmp(map, 't')

            % Loop over runs to load data
            counter = 0;
            con_dir = fullfile(pwd, '..', 'derivatives', subID, 'exp_glm1_norm_220/');

            for trial = 1:nTrials
                counter = counter + 1;

                % Load the contrast image for this run and trial
                con_file = fullfile(con_dir, sprintf('con_%04d.nii',trial));

                % Convert to CoSMoMVPA dataset
                ds = cosmo_fmri_dataset(con_file, ...
                    'mask', mask_fn, ... % Set brain mask
                    'targets', trial);

                % Store dataset
                datasets{counter} = ds;
            end

            % Combine all runs into a single dataset
            ds_per_run = cosmo_stack(datasets);

            % remove constant features
            ds_per_run = cosmo_remove_useless_data(ds_per_run);

            %% Make RDM

            % corrleate t maps
            rdm_corr = corr(ds_per_run.samples', 'type', 'Spearman');

        elseif strcmp(map, 'b')

            % get path
            betaPath = fullfile(pwd, '..', 'derivatives', subID, 'GLMsingleEstimates', ...
                'GLMsingle_betas.nii');

            % load beta map
            ds_per_run = cosmo_fmri_dataset(betaPath, ...
                'mask', mask_fn); % Set brain mask

            % add chunks and targets
            load(fullfile(pwd, '..', 'derivatives', subID, 'GLMsingleEstimates', ...
                'trialIDs.mat'));
            ds_per_run.sa.targets = trialIDs(:, 1);
            ds_per_run.sa.chunks = trialIDs(:, 2);

            % remove living room trials trials (target > 100)
            ds_per_run = cosmo_slice(ds_per_run, (ds_per_run.sa.targets <= 100), 1);

            % remove constant features
            ds=cosmo_remove_useless_data(ds_per_run);

            %% Get mean betas of glm single estimates

            meanBeta = [];
            for s1 = 1:nTrials

                % Subset data for the two stimuli
                ds_stim = cosmo_slice(ds, ds.sa.targets == s1);

                % take mean
                meanBeta(:, s1) = mean(ds_stim.samples)';

            end
            rdm_corr = corr(meanBeta, 'type', 'Spearman');


        else
            error('map not defined')
        end

        % Save the RDM
        subID2 = strrep(subID, '-', '');

        % get within and between category correlation
        rdm_corr(eye(nTrials, nTrials) == 1) = 0;
        withinCate = [squareform(rdm_corr(1:nTrials/2, 1:nTrials/2)), ...
            squareform(rdm_corr(nTrials/2 + 1:end, nTrials/2 + 1:end))];
        betweenCate = reshape(rdm_corr(nTrials/2 + 1:end, 1:nTrials/2), 1, []);

        % store in structure
        res.(map).(subID2).(mask_label_short).rdm_corr = rdm_corr;
        res.(map).(subID2).(mask_label_short).mean_cor = mean(squareform(rdm_corr));
        res.(map).(subID2).(mask_label_short).corr_cate_diff =  mean(withinCate) - mean(betweenCate);

    end
end

%% saving

% define output folder and name
outputFolder = fullfile(pwd, '..', 'derivatives', 'group_level', 'RDM');
outputName = ['results_RDM_of_mean_correaltion_on_', map, '-map.mat'];

% save 
save(fullfile(outputFolder, outputName), 'res')



%% ploting

% Get list of subjects and ROIs
subjects = fieldnames(res.(map));
masks = fieldnames(res.(map).(subjects{1}));

% Initialize data storage
num_subjects = numel(subjects);
num_rois = numel(masks);
all_data = nan(num_subjects, num_rois);
all_rdms = nan(num_subjects, num_rois, nTrials, nTrials);

% Collect data
for i_sub = 1:num_subjects
    subID = subjects{i_sub};
    for i_roi = 1:num_rois
        mask_label = masks{i_roi};
        all_data(i_sub, i_roi) = res.(map).(subID).(mask_label).corr_cate_diff;
        all_rdms(i_sub, i_roi, :, :) = res.(map).(subID).(mask_label).rdm_corr;
    end
end

% Compute mean and standard deviation for each ROI
mean_data = mean(all_data, 1);
std_data = std(all_data, 0, 1);
mean_rdm = squeeze(mean(all_rdms, 1));

%% Create bar plot with comparison of within and between category
% correlation
figure;
hold on;

% Bar plot with error bars
bar_handle = bar(mean_data, 'FaceColor', 'flat');
errorbar(1:num_rois, mean_data, std_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add horizontal line at chance level
yline(0, '--r', 'No Category difference', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');

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
ylabel('Correlation diff');
ylim([min(min(all_data))-0.05, max(max(all_data))+0.05])
title('Within - Between category pairwise correaltion');

hold off;

%% make tiled figure with mean rsm of each ROI
figure;
title('Mean pairwise correlation');
allRoiRDMs = nan(nTrials^2, num_rois);

for i_roi = 1:num_rois
    nexttile

    % get RDM for ROI
    roiRDM = squeeze(mean_rdm(i_roi, :, :));
    allRoiRDMs(:, i_roi) = reshape(roiRDM, [], 1);

    % plot RDM
    imagesc(roiRDM, [-1, 1])
    colorbar;
    title(masks{i_roi});
end

% get inter-roi correlation
nexttile
corrRois = corr(allRoiRDMs, 'type', 'Spearman');

imagesc(corrRois, [0.5, 1])
colorbar;
title('inter-ROI correlation');

% add ROI labels
xticks(1:num_rois);
xticklabels(masks);
yticks(1:num_rois);
yticklabels(masks);


end

function res = PairwiseROIsLDAclassifier(cfg)

% evaluate input
if ~isfield(cfg, 'func_roi_names'); cfg.func_roi_names = {'PPA', 'TOS', 'RSC', 'LOC'}; end

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

        % get current mask 
        mask_label=cfg.rois{j};

        disp(['Using mask ', mask_label]);
        disp(char(datetime))

        % get datasetn in Cosmo format
        ds = loadCosmoDataset(cfg, subID, mask_label);

        %% Pairwise decoding
        % Initialize RDM
        rdm = zeros(cfg.nTrials, cfg.nTrials);
        disp('')

        % define options
        opt.max_feature_count = 5000;

        if isempty(gcp('nocreate'))
            parpool(8);
        end
        nTrials = cfg.nTrials;
        parfor stim1 = 1:nTrials
            for stim2 = 1:nTrials
                if ~(stim2 > stim1)
                    continue
                end

                % Subset data for the two stimuli
                ds_stim = cosmo_slice(ds, ds.sa.targets == stim1 | ds.sa.targets == stim2);

                % Rename target
                ds_stim.sa.targets = (ds_stim.sa.targets == stim1) + 1;

                % Define partitions
                partitions=cosmo_nfold_partitioner(ds_stim);

                % get predictions for each fold
                [~ ,accuracy] = cosmo_crossvalidate(ds_stim, classifier, partitions, opt);

                % Store the accuracy
                rdm(stim2, stim1) = accuracy;
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
all_rdm_data = nan(cfg.nTrials, cfg.nTrials, num_subjects, num_rois);

% Collect data
for i_sub = 1:num_subjects
    subID = subjects{i_sub};
    for i_roi = 1:num_rois
        mask_label = masks{i_roi};
        all_data(i_sub, i_roi) = res.(subID).(mask_label).mean_accuracy;
        all_rdm_data(:, :, i_sub, i_roi) = res.(subID).(mask_label).rdm;
    end
end


%% bar plot
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
ylim([0.45,0.55])
title('Mean Pairwise Decoding Results Across ROIs');

hold off;

%% rdm

% take mean across subjects
mean_all_rdm_data = squeeze(mean(all_rdm_data, 3));


figure;
title('Mean pairwise decoding accuracy');

for i_roi = 1:num_rois
    nexttile

    % get RDM for ROI
    roiRDM = mean_all_rdm_data(:, :, i_roi);
    allRoiRDMs(:, i_roi) = reshape(roiRDM, [], 1);
    
    % plot RDM
    imagesc(roiRDM, [0.40, 0.60])
    colorbar;
    title(masks{i_roi});
end

% get inter-roi correlation
nexttile
corrRois = corr(allRoiRDMs, 'type', 'Spearman');

imagesc(corrRois, [-0.5, 0.5])
colorbar;
title('inter-ROI correlation');


% % get inter-roi correlation
% nexttile
% corrRois = corr(allRoiRDMs, 'type', 'Spearman');
%
% imagesc(corrRois, [0.5, 1])
% colorbar;
% title('inter-ROI correlation');
%
% % add ROI labels
% xticks(1:num_rois);
% xticklabels(masks);
% yticks(1:num_rois);
% yticklabels(masks);

%% Create bar plot with comparison of within and between category

diff_per_sub = nan(num_rois, num_subjects);
mean_diff = nan(1, num_rois);
std_diff = mean_diff;
for i_roi = 1:num_rois
    for i_sub = 1:num_subjects
        % get within and between category correlation
        roiRDM = all_rdm_data(:, :, i_sub, i_roi);
        withinCate = [squareform(roiRDM(1:cfg.nTrials/2, 1:cfg.nTrials/2)), ...
            squareform(roiRDM(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end))];
        betweenCate = reshape(roiRDM(cfg.nTrials/2 + 1:end, 1:cfg.nTrials/2), 1, []);
        diff_per_sub(i_roi, i_sub) = mean(withinCate) - mean(betweenCate);
    end
    % take the mean
    mean_diff(i_roi) = mean(diff_per_sub(i_roi, :));
    std_diff(i_roi) = std(diff_per_sub(i_roi, :));
end

% plot within - between difference
figure;
hold on;

% Bar plot with error bars
bar_handle = bar(mean_diff, 'FaceColor', 'flat');
errorbar(1:num_rois, mean_diff, std_diff, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add horizontal line at chance level
yline(0, '--r', 'No Category difference', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');

% Add jittered individual points
rng(0); % For reproducible jitter
jitter_amount = 0.1; % Adjust jitter spread
for i_roi = 1:num_rois
    x_jitter = i_roi + (rand(num_subjects, 1) - 0.5) * jitter_amount;
    scatter(x_jitter, diff_per_sub(i_roi, :), 30, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
    for i_sub = 1:num_subjects
        text(x_jitter(i_sub) + 0.04, diff_per_sub(i_roi, i_sub), subjects{i_sub}, 'FontSize', 8);
    end
end

% Customize plot
xticks(1:num_rois);
xticklabels(masks);
xlabel('ROI');
ylabel('Pariwise decoding diff');
ylim([min(min(diff_per_sub))-0.005, max(max(diff_per_sub))+0.005])
title('Within - Between category pairwise correaltion');

hold off;

%% is-rdm

figure;
title('IS RDMs across the whole RDM');

for i_roi = 1:num_rois
    nexttile

    % make a matrix with vectorized RDMs
    for i_sub = 1:num_subjects
        RDMmat(:, i_sub) = squareform(all_rdm_data(:, :, i_sub, i_roi));
    end

    % make and plot IS-RDM
    cfg.correlation_type = 'spearman';
    cfg.labels = cfg.subNums;
    cfg.cell_label_style = 'coef';
    cfg.plotting = true;
    cfg.new_figure = false;
    cfg.dissimilarity = false;
    [~, ~, ~] = make_RDM(RDMmat, cfg);
    title(masks{i_roi});
end


end
function res = CrossValdidatedRDM(cfg)

%% get data
dist = 'spearman';
% get number of ROI masks
nmasks=numel(cfg.rois);

% loop through subjects
for iSub = 1:length(cfg.subNums)

    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

    % progress report
    disp(['Starting cross validated RDMg for subject ',  num2str(cfg.subNums(iSub)), ' on ',...
        cfg.map, '-map']);

    for j=1:nmasks

        mask_label=cfg.rois{j};
        mask_label_short = split(mask_label, '.');
        mask_label_short = mask_label_short{1};

        disp(['Using mask ',  mask_label]);
        disp(char(datetime))

        % get datasetn in Cosmo format
        ds = loadCosmoDataset(cfg, subID, mask_label);

        % save data set
        save(fullfile(pwd, '..', 'derivatives', subID,...
            ['cosmo_dataset_', mask_label_short, '_', cfg.map, '_map_hpf_128.mat']), 'ds')

        %% Crossvalidated RDM

        % Initialize RDM
        rdm = zeros(cfg.nTrials, cfg.nTrials);

        % get combination to split 6 into 2 halfs
        splits = nchoosek(1:cfg.nRuns/2, cfg.nRuns/2/2);
        splits = splits(1:height(splits)/2, :);

        % just do even/odd for now
        %splits = [1,3,5,7,9];

        % abbaabbaab style
        %splits = [1,4,5,8,9];

        % equal number of even and odd
        %         for iSplit = 1:height(splits)
        %             current_split = splits(iSplit, :);
        %             numOdd = sum(mod(current_split, 2));
        %
        %             if numOdd < 2 || numOdd > 3
        %                 splits(iSplit, :) = nan;
        %             end
        %         end
        %         splits = splits(~isnan(splits(:,1)),:);

        % all splits rdm
        %all_splits_rdms = ones(cfg.nTrials, cfg.nTrials, height(splits));

        disp('')
        if isempty(gcp('nocreate'))
            parpool(8);
        end
        nTrials = cfg.nTrials;
        parfor stim1 = 1:cfg.nTrials
            disp(num2str(stim1))
            for stim2 = 1:nTrials
                if ~(stim2 > stim1)
                    continue
                end

                % Subset data for the two stimuli
                ds_stim = cosmo_slice(ds, ds.sa.targets == stim1 | ds.sa.targets == stim2);

                % Rename target
                ds_stim.sa.targets = (ds_stim.sa.targets == stim2) + 1;

                % loop over splits
                crossValdidatedRs = zeros(1, height(splits));
                for iSplit = 1:height(splits)
                    current_split = splits(iSplit, :);

                    % split data
                    dataSplit1Cond1 = cosmo_slice(ds_stim, ds_stim.sa.targets == 1 &...
                        ismember(ds_stim.sa.chunks, current_split));
                    dataSplit1Cond2 = cosmo_slice(ds_stim, ds_stim.sa.targets == 2 &...
                        ismember(ds_stim.sa.chunks, current_split));
                    dataSplit2Cond1 = cosmo_slice(ds_stim, ds_stim.sa.targets == 1 &...
                        ~ismember(ds_stim.sa.chunks, current_split));
                    dataSplit2Cond2 = cosmo_slice(ds_stim, ds_stim.sa.targets == 2 &...
                        ~ismember(ds_stim.sa.chunks, current_split));

                    % take mean for each split
                    meanData = zeros(width(dataSplit1Cond1.samples), 4);
                    meanData(:, 1) = mean(dataSplit1Cond1.samples);
                    meanData(:, 2) = mean(dataSplit1Cond2.samples);
                    meanData(:, 3) = mean(dataSplit2Cond1.samples);
                    meanData(:, 4) = mean(dataSplit2Cond2.samples);

                    % correlate them
                    distance = [];
                    if strcmp(dist, 'spearman')
                        rvals = corr(meanData, 'type', 'Spearman');
                        distance = rvals;
                    elseif strcmp(dist, 'euclidean')
                        distance = squareform(pdist(meanData'));
                    end

                    % get crossvalidation R
                    within = mean([distance(1, 3), distance(2, 4)]); % same cond, other split
                    between = mean([distance(1, 4), distance(2, 3)]); % other cond, other split
                    crossValdidatedRs(iSplit) =  within - between;

                    if isnan(crossValdidatedRs(iSplit))
                        warning('NaN has occured, dataset might be not suitable')
                    end
                end

                % Store the correlation
                rdm(stim2, stim1) = mean(crossValdidatedRs);

                % add to all msplits rdms
                %all_splits_rdms(stim2, stim1, :) = crossValdidatedRs;
            end
        end

        % make rdm symetric
        rdm = squareform(squareform(rdm));

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
all_rdm_data = nan(cfg.nTrials, cfg.nTrials, num_subjects, num_rois);

% Collect data
for i_sub = 1:num_subjects
    subID = subjects{i_sub};
    for i_roi = 1:num_rois
        mask_label = masks{i_roi};
        all_data(i_sub, i_roi) = res.(subID).(mask_label).mean_cor;
        all_rdm_data(:, :, i_sub, i_roi) = res.(subID).(mask_label).rdm;
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
yline(0, '--r', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');

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
ylabel('dissimilarity');
ylim([-0.02,0.02])
title('Mean of crossvalidated RDMs across ROIs');

hold off;

%% plot rdms 

% take mean across subjects
mean_all_rdm_data = squeeze(mean(all_rdm_data, 3));

figure;
title('Mean pairwise dissimilarity');

for i_roi = 1:num_rois
    nexttile

    % get RDM for ROI
    roiRDM = mean_all_rdm_data(:, :, i_roi);
    allRoiRDMs(:, i_roi) = reshape(roiRDM, [], 1);

    % plot RDM
    imagesc(roiRDM, [-0.05, 0.05])
    colorbar;
    title(masks{i_roi});
end

% get inter-roi correlation
nexttile
corrRois = corr(allRoiRDMs, 'type', 'Spearman');

imagesc(corrRois, [-0.5, 0.5])
colorbar;
title('inter-ROI correlation');


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

% plot within - between category difference
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
ylabel('Pariwise dissimilarity diff');
ylim([min(min(diff_per_sub))-0.0005, max(max(diff_per_sub))+0.0005])
title('Within - Between category pairwise crossvalidated dissimilarity');

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


% %% plot all RDMS for single sub and roi
%
% % Determine the number of splits
% nSplits = size(all_splits_rdms, 3);
%
% % Create a figure with tiled layout
% figure;
% %tiledlayout(7, 20, 'Padding', 'compact', 'TileSpacing', 'compact'); % Adjust grid size as needed
%
% for splitIdx = 1:nSplits
%     % Extract the RDM for the current split
%     rdm = all_splits_rdms(:, :, splitIdx);
%
%     % Convert the row of numbers in `splits` to a string
%     titleStr = mat2str(splits(splitIdx, :));
%
%     % Create a new tile for each RDM
%     nexttile;
%     imagesc(rdm, [-0.2, 0.4]); % Display the RDM with colorbar range
%
%     title(titleStr, 'FontSize', 8); % Add the split title
%     axis square; % Ensure square aspect ratio for RDM
% end
%
% % Add a global title or axis labels if needed
% sgtitle('RDMs for Each Split', 'FontSize', 14); % Super-title for the figure

end
% Configuration
cfg.n = length(cfg.subNums);  % Total number of subjects
cfg.nRuns = 10;      % Number of runs
categories = {'bathroom', 'kitchen'}; % Define categories
conditions = 4;

% Initialize results storage
ISC_diff = zeros(cfg.n, cfg.n, cfg.nRuns);

% loop through ROIs
for iRoi = 1:numel(cfg.rois)
    roi = char(cfg.rois{iRoi});
    mask_label_short = split(roi, '.');
    mask_label_short = mask_label_short{1};
    mask_label_short = mask_label_short(2:end);

    % Loop over subject pairs
    counter = 0;
    diffs = nan(1, nchoosek(cfg.n, 2));
    withinRDM = zeros(cfg.n, cfg.n);
    betweenRDM = zeros(cfg.n, cfg.n);
    for s1 = 1:cfg.n
        for s2 = s1:cfg.n
            if s1 == s2
                continue; % Skip same-subject comparisons
            end
            counter = counter +1;
            % Initialize run-wise correlation storage
            run_ISC_diff = zeros(1, cfg.nRuns);
            within = nan(1, cfg.n);
            between = nan(1, cfg.n);
            for iRun = 1:cfg.nRuns

                % Load time series for each subject and run
                for c = 1:conditions
                    if cfg.smoothing
                        fileName = ['smoothed_mean_timecourse_', categories{round(c/2)}, '_', ...
                            mask_label_short, '_run_', num2str(iRun), '.mat'];
                    else
                        fileName = ['mean_timecourse_', categories{round(c/2)}, '_', ...
                            mask_label_short, '_run_', num2str(iRun), '.mat'];
                    end

                    if mod(c, 2) == 1
                        subID = sprintf('sub-%03d',  cfg.subNums(s1));
                    else
                        subID = sprintf('sub-%03d',  cfg.subNums(s2));
                    end
                    % get timecourse
                    timecourseDir = fullfile(pwd, '..', 'derivatives', subID, 'timecourses');

                    load(fullfile(timecourseDir, fileName))
                    currentTimecourse = saveData;

                    if c == 1
                        pairData = zeros(length(saveData), 4);
                    end
                    pairData(:, c) = saveData';

                end

                % do correlation
                rvals = corr(pairData, 'type', 'Pearson', 'Rows', 'pairwise');
                within(iRun) = mean([rvals(1, 2), rvals(3, 4)]); % same cond, other split
                between(iRun) = mean([rvals(1, 4), rvals(2, 3)]); % other cond, other split

            end
            % write to RDM
            diffs(counter) =  mean(within, 'omitnan') - mean(between, 'omitnan');
            withinRDM(s1,s2) = mean(within, 'omitnan');
            withinRDM(s2,s1) = mean(within, 'omitnan');
            betweenRDM(s1,s2) = mean(between, 'omitnan');
            betweenRDM(s2,s1) = mean(between, 'omitnan');
        end
    end
    res.(mask_label_short).RDM = diffs;
    d.betweenRDM.(mask_label_short).RDM = betweenRDM;
    d.withinRDM.(mask_label_short).RDM = withinRDM;
    res.(mask_label_short).mean = mean(diffs, 'omitnan');
end

%% plotting
% Extract field names from the res structure (each corresponds to a mask)
mask_labels = fieldnames(res);

% Initialize arrays for storing means and individual data points
numMasks = numel(mask_labels);
means = zeros(1, numMasks);
allData = cell(1, numMasks);  % Store individual data points for each mask

% Extract data
for i = 1:numMasks
    mask = mask_labels{i};  % Get the current mask label
    means(i) = res.(mask).mean;  % Extract mean value
    allData{i} = res.(mask).RDM(:);  % Flatten RDM into a column vector
end

% Create bar plot
figure;
hold on;
bar(means, 'FaceColor', 'flat');  % Plot the mean values as bars
xticks(1:numMasks);
xticklabels(mask_labels);  % Set x-axis labels
ylabel('Mean Correlation');
title('ISC within - between categories');

% Scatter individual data points
for i = 1:numMasks
    x = i * ones(size(allData{i}));  % X-coordinates for scatter
    scatter(x, allData{i}, 50, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.1);  % Add jitter for visibility
end

% Improve appearance
grid on;
hold off;


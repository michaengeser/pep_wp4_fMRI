function d = neural_timecourse_intersub_cor(d, cfg)

warning('off');

% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'runSample'); cfg.runSample = 1:cfg.nRuns; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end
if ~isfield(cfg, 'detrend'); cfg.detrend = true; end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true; end
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = false; end
if ~isfield(cfg, 'smoothing'); cfg.smoothing = false; end

% convert to local configurations
if cfg.plotting; plotting = true; else plotting = false; end
if cfg.saving; saving = true; else saving = false; end
cfg.plotting = false;
cfg.saving = false;


% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % clear structure
    d.([category,'_RDM']).timecourseRDM = struct;

    % filter for category runs
    if strcmp(category, 'bathroom')
        runSample = cfg.runSample((mod(cfg.runSample,2) == 1));
    else
        runSample = cfg.runSample((mod(cfg.runSample,2) == 0));
    end

    % loop through ROIs
    for iRoi = 1:numel(cfg.rois)
        roi = char(cfg.rois{iRoi});
        mask_label_short = split(roi, '.');
        mask_label_short = mask_label_short{1};
        mask_label_short = mask_label_short(2:end);
        sub_table_concat = [];

        % loop through run
        runsCorrelationMatrix = nan(cfg.n, cfg.n, length(runSample));
        for iRun = 1:length(runSample)
            currentRun = runSample(iRun);

            % loop through subjects
            sub_table = [];
            for iSub = 1:cfg.n

                % get subjec id
                subID = sprintf('sub-%03d',  cfg.subNums(iSub));

                % get timecourse
                if cfg.cutTargets
                    timecourseDir = fullfile(cfg.outputPath, subID, 'timecourses');
                else
                    timecourseDir = fullfile(cfg.outputPath, subID, 'timecourses_with_targets');
                end

                if cfg.smoothing
                    fileName = ['smoothed_mean_timecourse_', category, '_', ...
                        mask_label_short, '_run_', num2str(currentRun), '.mat'];
                else
                    fileName = ['mean_timecourse_', category, '_', ...
                        mask_label_short, '_run_', num2str(currentRun), '.mat'];
                end

                load(fullfile(timecourseDir, fileName))
                currentTimecourse = saveData;

                if cfg.detrend
                    currentTimecourse = detrend(currentTimecourse);
                    currentTimecourse = currentTimecourse - mean(currentTimecourse);
                end

                % store in mat
                sub_table(:, iSub) = currentTimecourse';

            end

%             % get groub average
%             if cfg.regressOutMean
%                 % get mean
%                 groupMean = mean(sub_table, 2, 'omitnan');
% 
%                 % loop through subjects and regress out mean
%                 regressedTimecourses = zeros(size(sub_table)); % Initialize
%                 for iSub = 1:cfg.n
%                     % Design matrix: group-average timecourse and intercept
%                     X = [groupMean, ones(height(sub_table), 1)];
%                     % Perform regression
%                     beta = X \ sub_table(:, iSub); % Compute coefficients
%                     predicted = X * beta; % Predicted values based on the group average
%                     % Residual (subject timecourse with group average regressed out)
%                     regressedTimecourses(:, iSub) = sub_table(:, iSub) - predicted;
%                 end
% 
%                 % overwrite the subject table
%                 sub_table = regressedTimecourses;
%             end


            % make RDM
            cfg.labels = num2cell(cfg.subNums);
            cfg.correlation_type = 'Pearson';
            [~, runsCorrelationMatrix(:, :, iRun), ~] = make_RDM(sub_table, cfg);

            % store results
            % Add name, color and RDM fields to data struct
            d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).name = mask_label_short;
            d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).color = [0, 0, 0];
            d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).RDM = runsCorrelationMatrix(:, :, iRun);

            % concatanate tables across runs
            if isempty(sub_table_concat)
                sub_table_concat = sub_table;
            else
                sub_table_concat = [sub_table_concat; sub_table];
            end
        end


        % store results
        % Add name, color and RDM fields to data struct
        d.([category,'_RDM']).timecourseRDM(iRoi).name = mask_label_short;
        d.([category,'_RDM']).timecourseRDM(iRoi).color = [0, 0, 0];
        d.([category,'_RDM']).timecourseRDM(iRoi).RDM = mean(runsCorrelationMatrix, 3, 'omitnan');

        % make RDM for concatanated subtable
        cfg.labels = num2cell(cfg.subNums);
        cfg.correlation_type = 'Pearson';
        [~, concatMat, ~] = make_RDM(sub_table_concat, cfg);

        % store results
        % Add name, color and RDM fields to data struct
        d.([category,'_RDM']).timecourseRDM_concat(iRoi).name = mask_label_short;
        d.([category,'_RDM']).timecourseRDM_concat(iRoi).color = [0, 0, 0];
        d.([category,'_RDM']).timecourseRDM_concat(iRoi).RDM = concatMat;

    end % roi loop

    %% plotting
    if plotting

        % Create a new figure
        figure;
        tiledlayout(numel(cfg.rois), length(runSample), 'TileSpacing', 'Compact', 'Padding', 'Compact');

        % loop through rois and runs
        medianMat = nan(length(runSample), numel(cfg.rois));

        for iRoi = 1:numel(cfg.rois)
            for iRun = 1:length(runSample)

                % plot RDM for each subject and run
                nexttile;
                currentRDM = d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).RDM;
                if cfg.dissimilarity
                    imagesc(currentRDM, [0,2])
                else
                    imagesc(currentRDM, [-1,1])
                    colormap(cfg.colormaps.white_zero2)
                end
                xticklabels([]);
                yticklabels([]);

                % get median for each RDM
                currentRDM(eye(height(currentRDM)) == 1) = 0;
                medianValue = median(squareform(currentRDM), 'omitnan');
                if cfg.dissimilarity
                    medianMat(iRun, iRoi) = 1 - medianValue;
                end
            end
        end

        % add labels
        % Define row and column labels
        rowLabels = cfg.rois;
        colLabels = arrayfun(@num2str, 1:length(runSample), 'UniformOutput', false);

        % Add row labels (rois)
        for row = 1:numel(rowLabels)
            yPos = numel(rowLabels) + 1 - row - row / numel(rowLabels) / 2;
            xPos = -numel(colLabels) - 0.05;
            text(xPos, yPos, rowLabels{row},...
                'HorizontalAlignment', 'right','VerticalAlignment', 'middle',...
                'FontSize', 10, 'Units', 'normalized');
        end

        % Add column labels (numbers)
        for col = 1:numel(colLabels)
            yPos = -0.2;
            xPos = -numel(colLabels) + col - 0.3 + 0.08*col;
            text(xPos, yPos, colLabels{col},...
                'HorizontalAlignment', 'right','VerticalAlignment', 'middle',...
                'FontSize', 10, 'Units', 'normalized');
        end

        % add title
        sgtitle(['ISC on timecourse for ROIs and runs - ', category]);

        % bar plot
        figure;
        bar(mean(medianMat))
        ylabel('Median ISC')
        xlabel('ROI')
        xticklabels(cfg.rois)
        ylim([min(mean(medianMat))-0.01, max(mean(medianMat))+0.01])
        title(['Median ISC for each ROI (avergaed across runs) - ', category])

    end
end % category loop

warning('on');
end
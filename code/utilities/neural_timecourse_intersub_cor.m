function d = neural_timecourse_intersub_cor(d, cfg)

warning('off');

% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end

% convert to local configurations
if cfg.plotting; plotting = true; else plotting = false; end
if cfg.saving; saving = true; else saving = false; end
cfg.plotting = false;
cfg.saving = false;



% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % loop through ROIs
    for iRoi = 1:numel(cfg.rois)
        roi = char(cfg.rois{iRoi});
        mask_label_short = split(roi, '.');
        mask_label_short = mask_label_short{1};
        mask_label_short = mask_label_short(2:end);

        % loop through run
        runsCorrelationMatrix = nan(cfg.n, cfg.n, cfg.nRuns);
        for iRun = 1:cfg.nRuns

            % loop through subjects
            sub_table = [];
            for iSub = 1:cfg.n

                % get subjec id
                subID = sprintf('sub-%03d',  cfg.subNums(iSub));

                % get timecourse
                timecourseDir = fullfile(pwd, '..', 'derivatives', subID, 'timecourses');
                fileName = ['mean_timecourse_', category, '_', mask_label_short, ...
                    '_run_', num2str(iRun), '.mat'];
                load(fullfile(timecourseDir, fileName))
                currentTimecourse = saveData;

                % store in mat
                sub_table(:, iSub) = currentTimecourse';

            end

            % get goub average
            if cfg.regressOutMean
                % get mean
                groupMean = mean(sub_table, 2, 'omitnan');

                % loop through subjects and regress out mean
                regressedTimecourses = zeros(size(sub_table)); % Initialize
                for iSub = 1:cfg.n
                    % Design matrix: group-average timecourse and intercept
                    X = [groupMean, ones(height(sub_table), 1)];
                    % Perform regression
                    beta = X \ sub_table(:, iSub); % Compute coefficients
                    predicted = X * beta; % Predicted values based on the group average
                    % Residual (subject timecourse with group average regressed out)
                    regressedTimecourses(:, iSub) = sub_table(:, iSub) - predicted;
                end

                % overwrite the subject table
                sub_table = regressedTimecourses;
            end


            % make RDM
            cfg.labels = num2cell(cfg.subNums);
            cfg.correlation_type = 'Pearson';
            [~, runsCorrelationMatrix(:, :, iRun), ~] = make_RDM(sub_table, cfg);

            % store results
            % Add name, color and RDM fields to data struct
            d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).name = mask_label_short;
            d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).color = [0, 0, 0];
            d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).RDM = runsCorrelationMatrix(:, :, iRun);
        end


        % store results
        % Add name, color and RDM fields to data struct
        d.([category,'_RDM']).timecourseRDM(iRoi).name = mask_label_short;
        d.([category,'_RDM']).timecourseRDM(iRoi).color = [0, 0, 0];
        d.([category,'_RDM']).timecourseRDM(iRoi).RDM = mean(runsCorrelationMatrix, 3, 'omitnan');

    end % roi loop

    %% plotting
    if plotting

        % Create a new figure
        figure;
        tiledlayout(numel(cfg.rois), cfg.nRuns, 'TileSpacing', 'Compact', 'Padding', 'Compact');

        % loop through rois and runs
        medianMat = nan(cfg.nRuns, numel(cfg.rois));

        for iRoi = 1:numel(cfg.rois)
            for iRun = 1:cfg.nRuns

                % plot RDM for each subject and run
                nexttile;
                currentRDM = d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).RDM;
                imagesc(currentRDM, [-1,1])
                colormap(cfg.colormaps.white_zero)
                xticklabels([]);
                yticklabels([]);

                % get median for each RDM
                currentRDM(eye(height(currentRDM)) == 1) = 0;
                medianMat(iRun, iRoi) = median(squareform(currentRDM), 'omitnan'); 
            end
        end

        % add labels
        % Define row and column labels
        rowLabels = cfg.rois;
        colLabels = arrayfun(@num2str, 1:cfg.nRuns, 'UniformOutput', false);

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
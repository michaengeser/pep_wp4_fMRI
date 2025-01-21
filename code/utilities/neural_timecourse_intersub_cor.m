function d = neural_timecourse_intersub_cor(d, cfg, timecourses)

warning('off');

% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end

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

        % loop through run
        runsCorrelationMatrix = nan(cfg.n, cfg.n, cfg.nRuns);
        for iRun = 1:cfg.nRuns

            % loop through subjects
            sub_table = [];
            for iSub = 1:cfg.n

                % get subjec id
                subID = sprintf('sub%03d',  cfg.subNums(iSub));

                % get timecourse
                currentTimecourse = timecourses.(subID)(iRun).(mask_label_short).(['ROImeans', char(Category)]);

                % store in mat
                sub_table(:, iSub) = currentTimecourse';

            end

            % make RDM
            cfg.labels = num2cell(cfg.subNums);
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
        d.([category,'_RDM']).timecourseRDM(iRoi).RDM = mean(runsCorrelationMatrix, 3);


        if cfg.saving == 1
            % saving figure
            fig_path = fullfile(pwd, '..', 'figures', ['exp_', num2str(cfg.exp_num)], 'neural_ISC_RDMS', category);
            save_plot(fig_name, fig_path)
        end

    end % roi loop

    %% plotting
    if plotting

        % Create a new figure
        figure;
        t = tiledlayout(numel(cfg.rois), cfg.nRuns, 'TileSpacing', 'Compact', 'Padding', 'Compact');

        % loop through rois and runs
        medianMat = nan(cfg.nRuns, numel(cfg.rois));

        for iRoi = 1:numel(cfg.rois)
            for iRun = 1:cfg.nRuns

                % plot RDM for each subject and run
                nexttile;
                currentRDM = d.([category,'_RDM']).timecourseRDM(iRoi).runwiseRDM(iRun).RDM;
                imagesc(currentRDM, [-1,1])
                xticklabels([]);
                yticklabels([]);

                % get median for each RDM
                currentRDM(eye(height(currentRDM)) == 1) = 0;
                medianMat(iRun, iRoi) = median(squareform(currentRDM)); 
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
        ylim([0, max(mean(medianMat))+0.01])
        title(['Median ISC for each ROI (avergaed across runs) - ', category])

    end
end % category loop

warning('on');
end
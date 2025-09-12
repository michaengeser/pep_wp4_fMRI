function d = neural_timecourse_fc(d, cfg)

warning('off');

% defaults
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'saving'); cfg.saving = false; end
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'runSample'); cfg.runSample = 1:cfg.nRuns; end
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = false; end
if ~isfield(cfg, 'detrend'); cfg.detrend = true; end

plotting = cfg.plotting;
saving   = cfg.saving;

% loop through categories
for category = cfg.categories
    category = char(category);
    Category = strcat(upper(category(1)),lower(category(2:end))); % capitalize first letter

    % clear structure
    d.FC.(category).matrices = struct;

    % filter for category runs
    if strcmp(category, 'bathroom')
        runSample = cfg.runSample((mod(cfg.runSample,2) == 1));
    else
        runSample = cfg.runSample((mod(cfg.runSample,2) == 0));
    end

    % loop subjects
    for iSub = 1:cfg.n

        subID = sprintf('sub-%03d',  cfg.subNums(iSub));
        roi_timecourses = [];

        % collect all ROIs across runs
        for iRoi = 1:numel(cfg.rois)
            roi = char(cfg.rois{iRoi});
            mask_label_short = split(roi, '.');
            mask_label_short = mask_label_short{1};
            mask_label_short = mask_label_short(2:end);

            % concatenate across runs
            allRuns_tc = [];
            for iRun = 1:length(runSample)
                currentRun = runSample(iRun);

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
                allRuns_tc = [allRuns_tc, currentTimecourse];
            end

            % store for this ROI
            roi_timecourses(iRoi,:) = allRuns_tc;
        end

        % compute FC matrix for subject
        FCmat = corr(roi_timecourses');

        % store in struct
        d.FC.(category).matrices(iSub).name = subID;
        d.FC.(category).matrices(iSub).matrix = FCmat;
    end

    % compute mean across subjects
    allMats = cat(3, d.FC.(category).matrices.matrix);
    meanMat = mean(allMats,3, 'omitnan');

    d.FC.(category).meanMatrix = meanMat;

    %% plotting
    if plotting
        figure;
        imagesc(meanMat, [-1 1]);
        colorbar;
        title(['Mean functional connectivity - ', category]);
        xticks(1:numel(cfg.rois));
        yticks(1:numel(cfg.rois));
        xticklabels(cfg.rois);
        yticklabels(cfg.rois);
        xtickangle(45);
        colormap jet
    end

    %% Build Inter-Subject RDM from FC matrices

    % vectorize FC matrices (upper triangle only)
    mat_size = size(d.FC.(category).matrices(1).matrix);
    vecFC = nan(cfg.n, numel(squareform(ones(mat_size) - eye(mat_size))));
    for iSub = 1:cfg.n
        FCmat = d.FC.(category).matrices(iSub).matrix;
        vecFC(iSub,:) = 1-squareform(1-FCmat); % vectorize correaltion matrix
    end

    % subject-to-subject correlation
    fcCorrelationRDM = corr(vecFC', 'type', cfg.correlation_type, 'rows', 'pairwise');

    % convert correlation to dissimilarity RDM if desired
    if cfg.dissimilarity
        fcCorrelationRDM = 1 - fcCorrelationRDM;
    end

    % store results (one RDM for the whole ROI set)
    d.([category,'_RDM']).fcRDM(1).name  = 'allROIs';
    d.([category,'_RDM']).fcRDM(1).color = [0, 0, 0];
    d.([category,'_RDM']).fcRDM(1).RDM   = fcCorrelationRDM;
end

warning('on');
end

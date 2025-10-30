function d = compute_ROI_voxel_steps(cfg, d)

% defaults
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'runSample'); cfg.runSample = 1:cfg.nRuns; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end
if ~isfield(cfg, 'detrend'); cfg.detrend = true; end
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = false; end
if ~isfield(cfg, 'saving'); cfg.saving = true; end
if ~isfield(cfg, 'outputPath'); error('Need cfg.outputPath'); end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true; end
if ~isfield(cfg, 'runSample'); cfg.runSample = 1:cfg.nRuns; end
if ~isfield(cfg, 'numVoxels'); cfg.numVoxels = [1,2,5,10,20,50,100,200,500,1000,5000,inf]; end % inf = all voxels
if ~isfield(cfg, 'rois_of_interest'); cfg.rois_of_interest = {'LOC', 'PPA', 'TOS', 'LPFC'}; end

nSubs    = cfg.n;
subPairs = nchoosek(1:nSubs,2);
nPairs   = size(subPairs,1);

% loop ROIs
for iRoi = 1:numel(cfg.rois_of_interest)
    mask_label_short = cfg.rois_of_interest{iRoi};

    % loop categories
    for category = cfg.categories
        category = char(category);

        % check if file exists already
        outFile = fullfile(cfg.outputPath,'group_level', ...
            sprintf('ISC_%s_%s_voxel_steps.mat', mask_label_short, category));

        if cfg.skipIfExists && exist(outFile, 'file')
            disp(['Skipping existing file: ' outFile])
            continue
        end

        % select runs for this category
        if strcmp(category,'bathroom')
            runSample = cfg.runSample(mod(cfg.runSample,2)==1);
        else
            runSample = cfg.runSample(mod(cfg.runSample,2)==0);
        end

        % container for runwise results
        nROIs = length(cfg.rois_of_interest);
        d.voxel_steps_ISC.(category).(mask_label_short).runwiseVecRDMs = nan(nPairs, length(cfg.numVoxels), numel(runSample));

        % loop runs
        for iRun = 1:numel(runSample)
            currentRun = runSample(iRun);

            % preallocate subject × time × voxel steps
            sub_table = [];

            % load all subjects’ timecourses for this run
            for iSub = 1:nSubs
                subID = sprintf('sub-%03d',cfg.subNums(iSub));

                if cfg.cutTargets
                    tcDir = fullfile(cfg.outputPath,subID,'timecourses');
                else
                    tcDir = fullfile(cfg.outputPath,subID,'timecourses_with_targets');
                end

                fileName = fullfile(tcDir, ...
                    sprintf('mean_timecourses_voxel_steps_%s_%s_run-%02d.mat',mask_label_short, category, currentRun));

                tmp = load(fileName,'timecourses'); % [time × voxel steps]
                tc  = tmp.timecourses;

                if cfg.detrend
                    tc = detrend(tc);
                    tc = tc - mean(tc);
                end

                sub_table(:,:,iSub) = tc; % [time × steps × subject]
            end

            % loop through voxel number steps
            for iVox = 1:length(cfg.numVoxels)
                nVox = cfg.numVoxels(iVox);
                % subject × time for this voxel number
                tcMat = squeeze(sub_table(:,iVox,:)); % [time × subjects]

                % get groub average
                if cfg.regressOutMean
                    % get mean
                    groupMean = mean(tcMat, 2, 'omitnan');

                    % loop through subjects and regress out mean
                    regressedTimecourses = zeros(size(tcMat)); % Initialize
                    for iSub = 1:cfg.n
                        % Design matrix: group-average timecourse and intercept
                        X = [groupMean, ones(height(tcMat), 1)];
                        % Perform regression
                        beta = X \ tcMat(:, iSub); % Compute coefficients
                        predicted = X * beta; % Predicted values based on the group average
                        % Residual (subject timecourse with group average regressed out)
                        regressedTimecourses(:, iSub) = tcMat(:, iSub) - predicted;
                    end

                    % overwrite the subject table
                    tcMat = regressedTimecourses;
                end

                % correlation matrix across subjects
                R = corr(tcMat,'rows','pairwise');
                if cfg.dissimilarity
                    R = 1 - R;
                end

                % make everything NaN if too many NaNs
                if sum(isnan(R),'all') > sum(~isnan(R),'all')
                    R = nan(size(R));
                end
                R(eye(size(R) ) == 1) = 0;

                d.voxel_steps_ISC.(category).(mask_label_short).runwiseVecRDMs(:,iVox,iRun) = squareform(R)';

            end
        end

        % save runwise ISC results
        if cfg.saving

            meanISC = mean(d.voxel_steps_ISC.(category).(mask_label_short).runwiseVecRDMs, 3); % [pairs × steps × runs]
            d.voxel_steps_ISC.(category).(mask_label_short).mean = meanISC;

            save(outFile,'meanISC','subPairs','-v7.3');

        end
    end
end
end

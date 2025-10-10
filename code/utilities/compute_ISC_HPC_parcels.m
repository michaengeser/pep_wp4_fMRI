function d = compute_ISC_HPC_parcels(cfg, d)

% defaults
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'runSample'); cfg.runSample = 1:cfg.nRuns; end
if ~isfield(cfg, 'regressOutMean'); cfg.regressOutMean = true; end
if ~isfield(cfg, 'detrend'); cfg.detrend = true; end
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = false; end
if ~isfield(cfg, 'saving'); cfg.saving = true; end
if ~isfield(cfg, 'outputPath'); error('Need cfg.outputPath'); end
if ~isfield(cfg, 'TRstartBuffer'); cfg.TRstartBuffer = 3; end
if ~isfield(cfg, 'TRendBuffer'); cfg.TRendBuffer = 3; end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true; end
if ~isfield(cfg, 'searchlightSource'); cfg.searchlightSource = 'HPC'; end
if ~isfield(cfg, 'runSample'); cfg.runSample = 1:cfg.nRuns; end


nSubs    = cfg.n;
subPairs = nchoosek(1:nSubs,2);
nPairs   = size(subPairs,1);


% loop categories
for category = cfg.categories
    category = char(category);

    % select runs for this category
    if strcmp(category,'bathroom')
        runSample = cfg.runSample(mod(cfg.runSample,2)==1);
    else
        runSample = cfg.runSample(mod(cfg.runSample,2)==0);
    end

    % container for runwise results
    if strcmp(cfg.searchlightSource, 'HPC')
        nParcels = 180; % HPC atlas fixed
        d.HPC_ISC.(category).runwiseVecRDMs = nan(nPairs, nParcels, numel(runSample));
    elseif strcmp(cfg.searchlightSource, 'LPFC')
        if cfg.cutTargets
            tcDir = fullfile(cfg.outputPath, 'sub-101', 'timecourses');
        else
            tcDir = fullfile(cfg.outputPath, 'sub-101', 'timecourses_with_targets');
        end
        fileName = fullfile(tcDir, ...
            sprintf('mean_timecourses_voxel_steps_LPFC_%s_run-%02d.mat','bathroom',1));
        tmp = load(fileName,'step_sizes'); % [time × parcels]
        d.LPFC_steps.(category).runwiseVecRDMs = nan(nPairs, length(tmp.step_sizes), numel(runSample));
    end


    % loop runs
    for iRun = 1:numel(runSample)
        currentRun = runSample(iRun);

        % preallocate subject × time × parcel
        sub_table = [];

        % load all subjects’ timecourses for this run
        for iSub = 1:nSubs
            subID = sprintf('sub-%03d',cfg.subNums(iSub));

            if cfg.cutTargets
                tcDir = fullfile(cfg.outputPath,subID,'timecourses');
            else
                tcDir = fullfile(cfg.outputPath,subID,'timecourses_with_targets');
            end

            if strcmp(cfg.searchlightSource, 'HPC')
                fileName = fullfile(tcDir, ...
                    sprintf('mean_HPC_timecourses_%s_run-%02d.mat',category,currentRun));
            elseif strcmp(cfg.searchlightSource, 'LPFC')
                fileName = fullfile(tcDir, ...
                    sprintf('mean_timecourses_voxel_steps_LPFC_%s_run-%02d.mat',category,currentRun));
            end
            tmp = load(fileName,'timecourses'); % [time × parcels]
            tc  = tmp.timecourses;

            % remove TR buffers
            tc = tc(cfg.TRstartBuffer+1:end-cfg.TRendBuffer,:);

            if cfg.detrend
                tc = detrend(tc);
                tc = tc - mean(tc);
            end

            sub_table(:,:,iSub) = tc; % [time × parcels × subject]
        end

        % loop parcels
        for parc = 1:size(sub_table,2)
            % subject × time for this parcel
            tcMat = squeeze(sub_table(:,parc,:)); % [time × subjects]

            % regress out group mean if requested
            if cfg.regressOutMean
                groupMean = mean(tcMat,2,'omitnan');
                X = [groupMean, ones(size(groupMean))];
                for iSub = 1:nSubs
                    beta = X \ tcMat(:,iSub);
                    tcMat(:,iSub) = tcMat(:,iSub) - X*beta;
                end
            end

            % correlation matrix across subjects
            R = corr(tcMat,'rows','pairwise');
            if cfg.dissimilarity
                R = 1 - R;
            end

            % vectorize upper triangle

            if strcmp(cfg.searchlightSource, 'HPC')
                d.HPC_ISC.(category).runwiseVecRDMs(:,parc,iRun) = R(sub2ind(size(R),subPairs(:,1),subPairs(:,2)));
            elseif strcmp(cfg.searchlightSource, 'LPFC')
                d.LPFC_steps.(category).runwiseVecRDMs(:,parc,iRun) = R(sub2ind(size(R),subPairs(:,1),subPairs(:,2)));
            end

        end
    end

    % save runwise ISC results
    if cfg.saving

        if strcmp(cfg.searchlightSource, 'HPC')
            outFile = fullfile(cfg.outputPath,'group_level', ...
                sprintf('ISC_HPC_%s.mat',category));
            meanISC = mean(d.HPC_ISC.(category).runwiseVecRDMs, 3); % [pairs × parcels × runs]
            d.HPC_ISC.(category).mean = meanISC;

        elseif strcmp(cfg.searchlightSource, 'LPFC')
            outFile = fullfile(cfg.outputPath,'group_level', ...
                sprintf('ISC_LPFC_Steps_%s.mat',category));
            meanISC = mean(d.LPFC_steps.(category).runwiseVecRDMs, 3); % [pairs × parcels × runs]
            d.LPFC_steps.(category).mean = meanISC;

        end
        save(outFile,'meanISC','subPairs','-v7.3');

    end
end
end

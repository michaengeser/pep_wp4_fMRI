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

nSubs    = cfg.n;
subPairs = nchoosek(1:nSubs,2);
nPairs   = size(subPairs,1);
nParcels = 180; % HPC atlas fixed

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
    d.HPC_ISC.(category).runwiseVecRDMs = nan(nPairs, nParcels, numel(runSample));

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

            fileName = fullfile(tcDir, ...
                sprintf('mean_HPC_timecourses_%s_run-%02d.mat',category,currentRun));
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
        for parc = 1:nParcels
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
            d.HPC_ISC.(category).runwiseVecRDMs(:,parc,iRun) = R(sub2ind(size(R),subPairs(:,1),subPairs(:,2)));
        end
    end

    % save runwise ISC results
    if cfg.saving
        outFile = fullfile(cfg.outputPath,'group_level', ...
            sprintf('ISC_HPC_%s.mat',category));
        meanISC = mean(d.HPC_ISC.(category).runwiseVecRDMs, 3); % [pairs × parcels × runs]
        d.HPC_ISC.(category).mean = meanISC;
        save(outFile,'meanISC','subPairs','-v7.3');
    end
end
end

function d = compute_ISC_parcels(cfg, d)

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
if ~isfield(cfg, 'runSample'); cfg.runSample = 1:cfg.nRuns; end
if ~isfield(cfg, 'numVoxels'); cfg.numVoxels = [1,2,5,10,20,50,100,200,500,1000,5000,inf]; end % inf = all voxels
if ~isfield(cfg, 'rois_of_interest'); cfg.rois_of_interest = {'V1', 'LOC', 'PPA', 'TOS', 'LPFC'}; end

nSubs    = cfg.n;
subPairs = nchoosek(1:nSubs,2);
nPairs   = size(subPairs,1);
for iVox = 1:length(cfg.numVoxels)
    nVox = cfg.numVoxels(iVox);

    % loop categories
    for category = cfg.categories
        category = char(category);

        % check if file exists already

        outFile = fullfile(cfg.outputPath,'group_level', ...
            sprintf('ISC_HCP_%s_%d_voxels.mat', category, nVox));

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

        nParcels = 180; % HCP atlas fixed
        d.HCP_ISC.(category).(['vox',num2str(nVox)]).runwiseVecRDMs = nan(nPairs, nParcels, numel(runSample));

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
                    sprintf('mean_HCP_timecourses_%s_run-%02d_%d_voxels.mat',category,currentRun, nVox));

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

                % correlation matrix across subjects
                R = corr(tcMat,'rows','pairwise');
                if cfg.dissimilarity
                    R = 1 - R;
                end

                % make everything NaN if too many NaNs
                if sum(isnan(R),'all') > sum(~isnan(R),'all')
                    R = nan(size(R));
                end
                d.HCP_ISC.(category).(['vox',num2str(nVox)]).runwiseVecRDMs(:,parc,iRun) = squareform(R)';
            end
        end

        % save runwise ISC results
        if cfg.saving

            meanISC = mean(d.HCP_ISC.(category).(['vox',num2str(nVox)]).runwiseVecRDMs, 3); % [pairs × parcels × runs]
            d.HCP_ISC.(category).(['vox',num2str(nVox)]).mean = meanISC;

            save(outFile,'meanISC','subPairs','-v7.3');

        end
    end
end
end

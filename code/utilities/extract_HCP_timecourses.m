function extract_HCP_timecourses(cfg)

% Default settings
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'skipIfExists'); cfg.skipIfExists = false; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85; end
if ~isfield(cfg, 'nVols'); cfg.nVols = 188; end
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = false; end

subs = cfg.subNums;

% Load HCP atlas once
atlas = load_untouch_nii(fullfile(pwd, '..', 'MNI_ROIs', 'wHPC_atlas2.nii'));
atlas_img = atlas.img;

for iSub = 1:length(subs)

    disp(['Subject: ', num2str(subs(iSub))])

    subID = sprintf('sub-%0.3d', subs(iSub));

    % make make output folder doesn't exist yet
    if cfg.cutTargets
        outdir = fullfile(cfg.outputPath, subID, 'timecourses');
    else
        outdir = fullfile(cfg.outputPath, subID, 'timecourses_with_targets');
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    for iRun = 1:cfg.nRuns

        % get category of the run
        if mod(iRun, 2) == 1
            currentCat = 'bathroom';
        else
            currentCat = 'kitchen';
        end

        % Output file
        outFile = fullfile(outdir, ...
            sprintf('mean_HPC_timecourses_%s_run-%02d.mat', currentCat, iRun));
        if cfg.skipIfExists && exist(outFile, 'file')
            disp(['Skipping existing file: ' outFile])
            continue
        end

        % Locate functional file
        funcFile = fullfile(cfg.outputPath, subID, 'func', ...
            sprintf('wr%sxxxx_task-scenes_run-%d_bold.nii', subID, iRun));

        if ~exist(funcFile, 'file')
            error('Missing file: %s', funcFile)
        end

        % Load functional 4D data
        v = load_untouch_nii(funcFile);
        func = single(v.img); % X × Y × Z × T

        % Exclude buffer TRs
        includedTRs = true(1, size(func,4));
        includedTRs(1:cfg.TRstartBuffer) = false;
        includedTRs(end-cfg.TRendBuffer+1:end) = false;
        func = func(:,:,:,includedTRs);

        nTRs = size(func,4);
        nParcels = 180;
        timecourses = nan(nTRs, nParcels, 'single');

        % Efficient reshaping
        func2D = reshape(func, [], nTRs); % voxels × time

        for p = 1:nParcels
            voxIdx = atlas_img(:) == p;
            if any(voxIdx)
                roi_ts = mean(func2D(voxIdx,:), 1);
                timecourses(:,p) = roi_ts';
            end
        end

        save(outFile, 'timecourses', 'nTRs', 'nParcels')
        disp(['Saved: ' outFile])

    end % run
end % subject
end

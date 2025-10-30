function extract_parcel_timecourses(cfg)

% Default settings
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'skipIfExists'); cfg.skipIfExists = true; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85; end
if ~isfield(cfg, 'nVols'); cfg.nVols = 188; end
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = false; end
if ~isfield(cfg, 'numVoxels'); cfg.numVoxels = [1,2,5,10,20,50,100,200,500,1000,5000,inf]; end % inf = all voxels
if ~isfield(cfg, 'searchlightSource'); cfg.searchlightSource = 'HCP'; end

subs = cfg.subNums;
main_path = fullfile(pwd, '..', 'derivatives');
contrast_folder = 'loc_glm1_norm';


% Load atlas once
if strcmp(cfg.searchlightSource, 'HCP')
    atlas = load_untouch_nii(fullfile(pwd, '..', 'MNI_ROIs', 'wHCP_atlas2.nii'));
    atlas_img = atlas.img;
    cfg.stopTooLargeVoxelNum = false;
elseif strcmp(cfg.searchlightSource, 'allROI')
    atlas = load_untouch_nii(fullfile(pwd, '..', 'MNI_ROIs', 'allRoiMask.nii'));
    atlas_img = atlas.img;
    cfg.stopTooLargeVoxelNum = true;
end

for iVox = 1:length(cfg.numVoxels)
    nVox = cfg.numVoxels(iVox);

    for iSub = 1:length(subs)

        disp(['Subject: ', num2str(subs(iSub))])
        subID = sprintf('sub-%0.3d', subs(iSub));

        % Load contrast files
        if strcmp(cfg.searchlightSource, 'HCP')
            contrast_path = fullfile(main_path, subID, contrast_folder);
            contrast_nii_path = fullfile(contrast_path, 'con_0004.nii');
            if ~isfile(contrast_nii_path)
                warning(['Contrast file not found: ', contrast_nii_path]);
            end
            contrast_nii = load_untouch_nii(contrast_nii_path);
            contrast_data = contrast_nii.img;
            contrast_data = contrast_data(:);
        elseif strcmp(cfg.searchlightSource, 'allROI')
        end

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

            if strcmp(cfg.searchlightSource, 'HCP')
                outFile = fullfile(outdir, ...
                    sprintf('mean_HCP_timecourses_%s_run-%02d_%d_voxels.mat', currentCat, iRun, nVox));
            elseif strcmp(cfg.searchlightSource, 'allROI')
                outFile = fullfile(outdir, ...
                    sprintf('mean_allROI_timecourses_%s_run-%02d_%d_voxels.mat', currentCat, iRun, nVox));
            end

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
            %         includedTRs = true(1, size(func,4));
            %         includedTRs(1:cfg.TRstartBuffer) = false;
            %         includedTRs(end-cfg.TRendBuffer+1:end) = false;
            %         func = func(:,:,:,includedTRs);

            nTRs = size(func,4);
            nParcels = length(unique(atlas_img)) - 1; % all non-zero parcels
            timecourses = nan(nTRs, nParcels, 'single');

            % Efficient reshaping
            func2D = reshape(func, [], nTRs); % voxels × time

            for p = 1:nParcels
                voxIdx = atlas_img(:) == p;
                if any(voxIdx)

                    if nVox ~= inf

                        % Apply ROI mask to contrast map
                        masked_flat = contrast_data(voxIdx);
                        masked_flat(isnan(masked_flat)) = 0;

                        % Sort voxels by contrast
                        [tVals, sorted_indices] = sort(masked_flat,'descend');

                        if strcmp(cfg.searchlightSource, 'HCP')
                            sorted_indices(tVals <= 0) = []; % drop voxels wiht zero or negative t
                        end

                        if length(sorted_indices) >= nVox
                            top_idx = sorted_indices(1:nVox);
                        else
                            top_idx = sorted_indices;
                        end

                        voxIdx_values = find(voxIdx);
                        roi_ts = mean(func2D(voxIdx_values(top_idx),:), 1); % voxels x time
                        % end
                    else
                        roi_ts = mean(func2D(voxIdx,:), 1);
                    end
                    timecourses(:,p) = roi_ts';
                end
            end

            save(outFile, 'timecourses', 'nTRs', 'nParcels')
            disp(['Saved: ' outFile])

        end % run
    end % subject
end % n voxel loop
end

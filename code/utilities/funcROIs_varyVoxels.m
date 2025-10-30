function funcROIs_varyVoxels(cfg)
% Create functional ROIs with variable voxel counts and extract timecourses

%% Configuration
if ~isfield(cfg, 'func_roi_names'); cfg.func_roi_names = {'LPFC'}; end
if ~isfield(cfg, 'numVoxels'); cfg.numVoxels = [1,2,5,10,20,50,100,200,500,1000,5000,inf]; end % inf = all voxels
if ~isfield(cfg, 'fmri_pattern'); cfg.fmri_pattern = 'wsub-*_task-*_bold.nii'; end % adjust to your data
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'skipIfExists'); cfg.skipIfExists = false; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85; end
if ~isfield(cfg, 'nVols'); cfg.nVols = 188; end
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = false; end
if ~isfield(cfg, 'TRstartBuffer'); cfg.TRstartBuffer = 3; end
if ~isfield(cfg, 'TRendBuffer'); cfg.TRendBuffer = 3; end

% Paths
main_path = fullfile(pwd, '..', 'derivatives');
contrast_folder = 'loc_glm1_norm';
contrast_files = {'con_0001.nii', 'con_0003.nii', 'con_0004.nii'};
ROI_dir = fullfile(pwd, '..', 'MNI_ROIs');
output_dir = fullfile(ROI_dir, 'func_ROIs');

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% Load ROI masks
roi_masks = cell(1, numel(cfg.func_roi_names));
for j = 1:numel(cfg.func_roi_names)
    roiFilName = ['w', cfg.func_roi_names{j}, '.nii'];
    mask_nii = load_untouch_nii(fullfile(ROI_dir, roiFilName));
    newMaskImg = mask_nii.img;
    if max(newMaskImg(:)) > 1
        newMaskImg = newMaskImg / max(newMaskImg(:));
    end
    roi_masks{j} = newMaskImg;
end


% Generate ROIs for each region
for iRoi = 1:numel(cfg.func_roi_names)

    %% Process each subject
    for iSub = 1:numel(cfg.subNums)
        subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
        fprintf('Processing subject %s...\n', subID);

        % make make output folder doesn't exist yet
        if cfg.cutTargets
            outdir = fullfile(cfg.outputPath, subID, 'timecourses');
        else
            outdir = fullfile(cfg.outputPath, subID, 'timecourses_with_targets');
        end
        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end

        % Load contrast maps
        contrast_path = fullfile(main_path, subID, contrast_folder);

        contrast_data = cell(1, numel(contrast_files));
        for c = 1:numel(contrast_files)
            contrast_nii_path = fullfile(contrast_path, contrast_files{c});
            if ~isfile(contrast_nii_path)
                warning(['Contrast file not found: ', contrast_nii_path]);
                continue;
            end
            contrast_nii = load_untouch_nii(contrast_nii_path);
            contrast_data{c} = contrast_nii.img;
        end

        % Pick contrast based on ROI
        if ismember(cfg.func_roi_names{iRoi}, {'PPA','TOS','RSC'})
            contrast_map = single(contrast_data{1});
        elseif ismember(cfg.func_roi_names{iRoi}, {'LOC'})
            contrast_map = single(contrast_data{2});
        elseif ismember(cfg.func_roi_names{iRoi}, {'LPFC'})
            contrast_map = single(contrast_data{3});
        else
            continue;
        end

        roi_mask = roi_masks{iRoi};

        % Apply ROI mask to contrast map
        masked_contrast = single(contrast_map) .* single(roi_mask);
        masked_flat = masked_contrast(:);
        masked_flat(isnan(masked_flat)) = 0;

        % Sort voxels by contrast
        [~, sorted_indices] = sort(masked_flat,'descend');
        sorted_indices(masked_flat(sorted_indices)==0) = []; % drop non-ROI voxels

        for iRun = 1:cfg.nRuns

            % get category of the run
            if mod(iRun, 2) == 1
                currentCat = 'bathroom';
            else
                currentCat = 'kitchen';
            end

            % Output file
            outFile = fullfile(outdir, ...
                sprintf('mean_timecourses_voxel_steps_%s_%s_run-%02d.mat', ...
                cfg.func_roi_names{iRoi}, currentCat, iRun));
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

            % reshpae func
            nTRs = size(func,4);
            func_2D = reshape(func, [], nTRs);

            % Preallocate matrix (time x voxels) for each step
            nSizes = length(cfg.numVoxels);
            timecourses = nan(nTRs, nSizes, 'single');
            for iStep = 1:nSizes
                step = cfg.numVoxels(iStep);

                % if step if larger than ROI size
                if length(sorted_indices) <= step
                    top_idx = sorted_indices;
                else
                    top_idx = sorted_indices(1:step);
                end
                roi_tc = func_2D(top_idx,:); % voxels x time
                if step > 1
                    timecourses(:, iStep) = mean(roi_tc)';  % take mean
                elseif step == 1
                    timecourses(:, iStep) = roi_tc';  % take mean
                end
            end

            save(outFile, 'timecourses', 'nTRs')
            disp(['Saved: ' outFile])

        end
    end
end

disp('All subjects processed.');
end

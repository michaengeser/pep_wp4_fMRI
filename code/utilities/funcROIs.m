function funcROIs(cfg)

%% Configuration
if ~isfield(cfg, 'func_roi_names'); cfg.func_roi_names = {'PPA', 'TOS', 'RSC', 'LOC', 'LPFC'}; end
if ~isfield(cfg, 'n_voxels'); cfg.n_voxels = 200; end

% Define subject IDs
main_path = fullfile(pwd, '..', 'derivatives'); % Base directory for derivatives
contrast_folder = 'loc_glm1_norm'; % Subfolder containing contrasts
contrast_files = {'con_0001.nii', 'con_0003.nii', 'con_0004.nii'}; % Scene > Objects, Objects > Scramble, , Intact > Scramble
ROI_dir = fullfile(pwd, '..', 'MNI_ROIs'); % Directory where anatomical ROI are
output_dir = fullfile(ROI_dir, 'func_ROIs'); % Directory to save functional ROIs

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Load ROI masks
roi_masks = cell(1, numel(cfg.func_roi_names));
for j = 1:numel(cfg.func_roi_names)
    roiFilName = ['w', cfg.func_roi_names{j}, '.nii'];
    mask_nii = load_untouch_nii(fullfile(ROI_dir, roiFilName));
    newMaskImg = mask_nii.img;
    if max(max(max(double(newMaskImg)))) > 1
        newMaskImg = newMaskImg/max(max(max(double(newMaskImg))));
    end
    roi_masks{j} = newMaskImg;
end

%% Process each subject
for s = 1:numel(cfg.subNums)
    sub_id = sprintf('sub-%0.3d', cfg.subNums(s));
    %disp(['Processing subject: ', sub_id]);

    % Subject-specific contrast path
    contrast_path = fullfile(main_path, sub_id, contrast_folder);

    % Load contrast maps
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

    % Generate ROIs for each region
    for r = 1:numel(cfg.func_roi_names)
        %disp(['  Generating ROI: ', cfg.func_roi_names{r}]);

        % Select the appropriate contrast map based on ROI
        if ismember(cfg.func_roi_names{r}, {'PPA', 'TOS', 'RSC'})
            contrast_map = single(contrast_data{1}); % Scene > Objects
        elseif ismember(cfg.func_roi_names{r}, {'LOC'})
            contrast_map = single(contrast_data{2}); % Objects > Scramble
        elseif ismember(cfg.func_roi_names{r}, {'LPFC'})
            contrast_map = single(contrast_data{3}); % Intact > Scramble
        else
            continue;
        end

        % Get the corresponding ROI mask
        roi_mask = roi_masks{r};

        % Apply the mask to the contrast map
        masked_contrast = single(contrast_map) .* single(roi_mask);

        % Flatten the masked contrast map for voxel selection
        masked_flat = masked_contrast(:);
        masked_flat(isnan(masked_flat)) = 0;

        % Sort voxels by intensity
        [sorted_values, sorted_indices] = sort(masked_flat, 'descend');

        % Select the top N voxels
        top_voxel_indices = sorted_indices(1:cfg.n_voxels);
        functional_roi = zeros(size(masked_contrast));
        functional_roi(top_voxel_indices) = 1;

        % Save the functional ROI as a new NIfTI file
        sub_output_dir = fullfile(output_dir, sub_id);
        if ~exist(sub_output_dir, 'dir')
            mkdir(sub_output_dir);
        end
        output_fn = fullfile(sub_output_dir, [cfg.func_roi_names{r}, '_funcROI.nii']);
        roi_nii = mask_nii; % Use the mask NIfTI as a template
        roi_nii.img = functional_roi;
        save_untouch_nii(roi_nii, output_fn);

        disp(['    Saved functional ROI for ', cfg.func_roi_names{r}, ' to ', output_fn]);
    end
end

disp('All subjects processed.');
end
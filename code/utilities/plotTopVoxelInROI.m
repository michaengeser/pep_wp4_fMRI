
%% Configuration
if ~isfield(cfg, 'roi_names'); cfg.roi_names = {'PPA', 'TOS', 'RSC', 'LOC'}; end
if ~isfield(cfg, 'n_voxels'); cfg.n_voxels = 100; end

% Define subject IDs
main_path = fullfile(pwd, '..', 'derivatives'); % Base directory for derivatives
contrast_folder = 'loc_glm1_norm'; % Subfolder containing contrasts
contrast_files = {'con_0001.nii', 'con_0003.nii'}; % Scene > Objects, Objects > Scramble
ROI_dir = fullfile(pwd, '..', 'MNI_ROIs'); % Directory where anatomical ROI are
output_dir = fullfile(ROI_dir, 'func_ROIs'); % Directory to save functional ROIs

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Load ROI masks
roi_masks = cell(1, numel(cfg.roi_names));
for j = 1:numel(cfg.roi_names)
    roiFilName = ['w', cfg.roi_names{j}, '.nii'];
    mask_nii = load_untouch_nii(fullfile(ROI_dir, roiFilName));
    newMaskImg = mask_nii.img;
    if max(max(max(double(newMaskImg)))) > 1
        newMaskImg = newMaskImg/max(max(max(double(newMaskImg))));
    end
    roi_masks{j} = newMaskImg;
end

%% Initialize structure to store top 100 voxel values
top_voxel_values = struct();
for r = 1:numel(cfg.roi_names)
    top_voxel_values.(cfg.roi_names{r}) = [];
end

%% Process each subject
for s = 1:numel(cfg.subNums)
    sub_id = sprintf('sub-%0.3d', cfg.subNums(s));
    disp(['Processing subject: ', sub_id]);

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
    for r = 1:numel(cfg.roi_names)
        disp(['  Generating ROI: ', cfg.roi_names{r}]);

        % Select the appropriate contrast map based on ROI
        if ismember(cfg.roi_names{r}, {'PPA', 'TOS', 'RSC'})
            contrast_map = single(contrast_data{1}); % Scene > Objects
        elseif ismember(cfg.roi_names{r}, {'LOC'})
            contrast_map = single(contrast_data{2}); % Objects > Scramble
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

        % Save the top 100 voxel values for plotting
        top_100_values = sorted_values(1:100);
        top_voxel_values.(cfg.roi_names{r}) = [top_voxel_values.(cfg.roi_names{r}); top_100_values];

        % Save the functional ROI as a new NIfTI file
        sub_output_dir = fullfile(output_dir, sub_id);
        if ~exist(sub_output_dir, 'dir')
            mkdir(sub_output_dir);
        end
        output_fn = fullfile(sub_output_dir, [cfg.roi_names{r}, '_funcROI.nii']);
        roi_nii = mask_nii; % Use the mask NIfTI as a template
        roi_nii.img = functional_roi;
%         save_untouch_nii(roi_nii, output_fn);

        disp(['    Saved functional ROI for ', cfg.roi_names{r}, ' to ', output_fn]);
    end
end

%% Generate scatter plots for top 100 voxel values
fig_dir = fullfile(output_dir, 'plots');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

for r = 1:numel(cfg.roi_names)
    roi_name = cfg.roi_names{r};
    values = top_voxel_values.(roi_name);
    
    if isempty(values)
        warning(['No data available for ROI: ', roi_name]);
        continue;
    end

    % Prepare x-axis values (subjects)
    x_values = repmat(1:numel(cfg.subNums), 100, 1);
    x_values = x_values(:);
    
    % Prepare y-axis values (t-values)
    y_values = values(:);

    % Create scatter plot
    figure;
    scatter(x_values, y_values, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('Subjects');
    ylabel('Voxel t-values');
    title(['Top 100 voxel t-values for ', roi_name]);
    grid on;

    % Save the figure
%     fig_path = fullfile(fig_dir, [roi_name, '_scatter_plot.png']);
%     saveas(gcf, fig_path);
%     close(gcf);

%     disp(['Saved scatter plot for ', roi_name, ' at ', fig_path]);
end

disp('All subjects processed and plots generated.');


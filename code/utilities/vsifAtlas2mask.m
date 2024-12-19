% Define paths
atlas_path = fullfile(pwd, '..', 'MNI_ROIs', 'visfAtlas_MNI152_volume.nii');
output_dir = fullfile(pwd, '..', 'MNI_ROIs');

% Load the atlas
atlas = niftiread(atlas_path);
atlas_info = niftiinfo(atlas_path);

% Define ROI indices for each mask


% roi_indices = struct(...
%     'lV1', [12, 15], ... % lh_v1d_retinotopic (12) and lh_v1v_retinotopic (15)
%     'rV1', [28, 31], ... % rh_v1d_retinotopic (28) and rh_v1v_retinotopic (31)
%     'lV2', [13, 16], ... % lh_v2d_retinotopic (13) and lh_v2v_retinotopic (16)
%     'rV2', [29, 32]);    % rh_v2d_retinotopic (29) and rh_v2v_retinotopic (32)

roi_indices = struct(...
    'V1', [12, 15, 28, 31], ...
    'V2', [13, 16, 29, 32]);

% roi_indices = struct(...
%     'visualCortex', 1:33);   % get all visual ROIs

% Loop through ROIs and create masks
fields = fieldnames(roi_indices);
for i = 1:numel(fields)
    roi_name = fields{i};
    roi_values = roi_indices.(roi_name);

    % Create binary mask
    binary_mask = ismember(atlas, roi_values);

    % Modify info
    new_info = atlas_info;
    mask_file_name =  [roi_name '.nii'];
    new_info.Filename = strrep(atlas_info.Filename,...
        'visfAtlas_MNI152_volume.nii', mask_file_name);
    new_info.Filemoddate = char(datetime);

    % Save the binary mask
    output_path = fullfile(output_dir, mask_file_name);
    niftiwrite(int16(binary_mask), output_path, new_info);

    fprintf('Saved %s at %s\n', roi_name, output_path);
end

disp('All ROI masks have been created and saved.');

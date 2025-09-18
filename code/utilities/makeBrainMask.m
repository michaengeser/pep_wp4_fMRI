function makeBrainMask(cfg)

if ~isfield(cfg, 'threshold'); cfg.threshold = 0.5; end

% Loop through remaining subjects
for iSub = 1:cfg.n

    % get subject pmask
    subID = sprintf('sub-%03d',  cfg.subNums(iSub));
    subject_mask_path = fullfile(pwd, '..', 'derivatives', subID, 'loc_glm1_norm', 'mask.nii');
    try
        nii = load_untouch_nii(subject_mask_path);
    catch
        continue
    end
    subject_mask = double(nii.img);

    % Create intersection: keep only voxels that are 1 in all subjects
    if iSub == 1
        group_intersection_mask = double(nii.img); % Load the first subject's mask as reference
        mask_sum = zeros(size(nii.img)); % Initialize sum mask
    else
        group_intersection_mask = group_intersection_mask & subject_mask;
        mask_sum = mask_sum + double(nii.img);
    end
end

% Convert back to uint8 (binary mask)
group_intersection_mask = uint8(group_intersection_mask);

% Save the final group mask
nii.img = group_intersection_mask;
save_path = fullfile(pwd, '..', 'derivatives', 'group_level');
if ~exist(save_path, 'dir')
    mkdir(save_path);
end
save_untouch_nii(nii, fullfile(save_path, 'group_intersection_mask.nii'));

% Apply threshold
voxel_threshold = round(cfg.threshold * cfg.n);

% Create final mask
group_tresholded_mask = uint8(mask_sum >= voxel_threshold);

% Save final mask
nii.img = group_tresholded_mask;
save_untouch_nii(nii, fullfile(save_path, 'group_mask_thresholded.nii'));

% Create final mask
group_mask_common = uint8(mask_sum >= 1);

% Save final mask
nii.img = group_mask_common;
save_untouch_nii(nii, fullfile(save_path, 'group_mask_common.nii'));


disp('Group whole-brain mask saved as group_mask.nii');

end

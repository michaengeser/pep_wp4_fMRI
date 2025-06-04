function cfg = viewMasks(cfg)

sum_mask = [];
for iSub = 1:cfg.n
    subID = sprintf('sub-%03d',  cfg.subNums(iSub));
    filePath = (fullfile(pwd, '..', 'MNI_ROIs', 'func_ROIs', subID, ...
        'LPFC_funcROI.nii'));
    nii = load_untouch_nii(filePath);
    if isempty(sum_mask)
        sum_mask = double(nii.img > 0);
    else
        sum_mask = sum_mask + double(nii.img > 0); % binary mask
    end
end

output_nii = nii; % use one subject's header
output_nii.img = sum_mask;
output_nii.hdr.dime.datatype = 16; % float32
output_nii.hdr.dime.bitpix = 32;
output_path = (fullfile(pwd, '..', 'MNI_ROIs', 'func_ROIs', ...
        'group_pfc_visual_response_map.nii'));
save_untouch_nii(output_nii, output_path);

end
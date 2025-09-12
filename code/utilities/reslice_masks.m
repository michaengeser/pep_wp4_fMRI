clc
clear

%% create a functional image for the space to reslice to
matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(pwd, '..', 'derivatives', 'sub-101', 'loc_glm1_norm', 'mask.nii')};

%% bring MNI mask into voxel space
matlabbatch{1}.spm.spatial.coreg.write.source = {
    fullfile(pwd, '..', 'MNI_ROIs', 'lPPA.img,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'lTOS.img,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'lRSC.img,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'lLOC.img,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'lV1.nii,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'lV2.nii,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'V2.nii,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'V1.nii,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'rV1.nii,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'rV2.nii,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'rTOS.img,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'rPPA.img,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'rRSC.img,1')
    fullfile(pwd, '..', 'MNI_ROIs', 'rLOC.img,1')
    };
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4; %%% MAYBE CHANGE TO 0
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'w';
spm_jobman('run_nogui',matlabbatch)
clear matlabbatch




%% create a functional image for the space to reslice to
matlabbatch{1}.spm.spatial.coreg.write.ref = {fullfile(pwd, '..', 'derivatives', 'sub-101', 'loc_glm1_norm', 'mask.nii')};

%% bring MNI mask into voxel space
matlabbatch{1}.spm.spatial.coreg.write.source = {fullfile(pwd, '..', 'MNI_ROIs', 'aal.nii,1')};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4; %%% MAYBE CHANGE TO 0
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'w';
spm_jobman('run_nogui',matlabbatch)
clear matlabbatch


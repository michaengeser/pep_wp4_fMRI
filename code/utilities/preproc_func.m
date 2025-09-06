function preproc_func(cfg)


%% Initialize SPM
matlabbatch = [];
spm('defaults', 'fmri');
spm_jobman('initcfg');

%% Run preprocessing
for iSub = 1:length(cfg.subNums)

    % subject ID
    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
    subIDx = sprintf('sub-%0.3dxxxx', cfg.subNums(iSub));

    % check if preprocessed files exit already
    outputDir = fullfile(cfg.outputPath, subID, 'func');
    if exist(outputDir, 'dir')
        motionResgressors = fullfile(outputDir, sprintf('rp_%sxxxx_task-scenes_run-%d_bold_%0.5d.txt', subID, cfg.nRuns, 1));
        if exist(motionResgressors, 'file')
            disp(['Preprocessed data for subject ', num2str(cfg.subNums(iSub)), ' exists already'])
            continue
        end 
    else
        mkdir(outputDir);
    end

    % define path to T1
    anatomy = fullfile(cfg.sourcedataPath, subID, 'anat', sprintf('%s_T1w.nii', subIDx));

    % define paths to BOLD files
    funcFiles = {};
    for iRun = 1:cfg.nRuns + 1
        if iRun == 1

            % get localizer
            singleRunBOLD = fullfile(cfg.sourcedataPath, subID, 'func', ...
                sprintf('%s_task-localizer_bold.nii', subIDx));
        else
            % get experimental runs
            singleRunBOLD = fullfile(cfg.sourcedataPath, subID, 'func', ...
                sprintf('%s_task-scenes_run-%d_bold.nii', subIDx, iRun-1));
        end

        % put localizer and experimental runs in unified cell array
        if exist('singleRunBOLD', 'var')
            funcFiles{iRun, 1} = singleRunBOLD;
        else
            warning(['File ', sprintf('%s_task-scenes_run-%d_bold.nii', subIDx, iRun-1),...
                ' is missing'])
        end
    end

    %% Define pre-processing steps
    for i = 1:length(funcFiles)
        % 4D to 3D-file conversion
        matlabbatch{i}.spm.util.split.vol = {strcat(funcFiles{i}, ',1')};
        matlabbatch{i}.spm.util.split.outdir = {outputDir};
    end

    % realignment
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{3}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{4}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{5}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{6}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{7}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{8}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{9}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{10}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{11}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{12}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.data{13}(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{13}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{i+1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{i+1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{i+1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{i+1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{i+1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{i+1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{i+1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{i+1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{i+1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{i+1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{i+1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{i+1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

    % coregistration
    matlabbatch{i+2}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{i+2}.spm.spatial.coreg.estimate.source = {strcat(anatomy, ',1')};
    matlabbatch{i+2}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{i+2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{i+2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{i+2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{i+2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % grey/white-matter segmentation
    matlabbatch{i+3}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{15}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{i+3}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{i+3}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{i+3}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(1).tpm = {fullfile(cfg.spmPath,'TPM.nii,1')};
    matlabbatch{i+3}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{i+3}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(2).tpm = {fullfile(cfg.spmPath,'TPM.nii,2')};
    matlabbatch{i+3}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{i+3}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(3).tpm = {fullfile(cfg.spmPath,'TPM.nii,3')};
    matlabbatch{i+3}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{i+3}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(4).tpm = {fullfile(cfg.spmPath,'TPM.nii,4')};
    matlabbatch{i+3}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{i+3}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(5).tpm = {fullfile(cfg.spmPath,'TPM.nii,5')};
    matlabbatch{i+3}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{i+3}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(6).tpm = {fullfile(cfg.spmPath,'TPM.nii,6')};
    matlabbatch{i+3}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{i+3}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{i+3}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{i+3}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{i+3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{i+3}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{i+3}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{i+3}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{i+3}.spm.spatial.preproc.warp.write = [1 1];
    matlabbatch{i+3}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{i+3}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
        NaN NaN NaN];

    % normalize functional runs
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(3) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 3)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(4) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 4)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(5) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 5)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(6) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 6)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(7) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 7)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{7}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(8) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 8)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{8}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(9) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 9)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{9}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(10) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 10)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{10}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(11) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 11)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{11}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(12) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 12)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{12}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.subj.resample(13) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 13)', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{13}, '.','rfiles'));
    matlabbatch{i+4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{i+4}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{i+4}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{i+4}.spm.spatial.normalise.write.woptions.prefix = 'w';

    % normalize anatomy
    matlabbatch{i+5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{i+5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{15}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{i+5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{i+5}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{i+5}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{i+5}.spm.spatial.normalise.write.woptions.prefix = 'w';

    %% Run pre-processing
    spm_jobman('run', matlabbatch)
    clear matlabbatch
    disp(['Done pre-processing subject ', char(cfg.subNums(iSub)), '!']);
end
end

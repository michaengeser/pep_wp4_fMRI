function level2_glm_func(subs)    
    %% Initialize SPM
    matlabbatch = [];
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    
    %% Define subjects and paths
    mainPath = '/Users/philippflieger/PhD/projects/fusion-study_scene-attractiveness/';
    fmriPath = fullfile(mainPath, 'sourcedata/fmri/bids');
    outPath = fullfile(fmriPath, 'derivatives/level2_glm');
    
    %% 2nd-level GLM: Localizer
    numCons = 3;
    
    %{
        --> Why 3 contrasts?
        1. scenes > objects
        2. scenes > scramble
        3. objects > scramble
    %}
    
    for iCon = 1:numCons
        conFolderName = ['con' num2str(iCon)];
    
        % create folder for contrast if it doesn't exist
        if exist(fullfile(outPath, 'loc', conFolderName)) > 0
            rmdir(fullfile(outPath, 'loc', conFolderName), 's');
        end
    
        mkdir(fullfile(outPath, 'loc', conFolderName));
    
        % specify model
        matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(outPath, 'loc', conFolderName)};
        conFiles = cell(length(subs), 1);
    
        % get each subject's files for the contrast of interest
        for iSub = 1:subs
            % get subject ID
            if length(num2str(iSub)) < 2
                subID = ['sub-00' num2str(iSub) 'xxxx'];
            else
                subID = ['sub-0' num2str(iSub) 'xxxx'];
            end
    
            % find contrast file
            conFiles{iSub,1} = fullfile( ...
                fmriPath, ...
                subID, ...
                sprintf( ...
                    'derivatives/loc_glm1_norm/con_000%s.nii,1', num2str(iCon) ...
                ) ...
            );
        end
    
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = conFiles;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
        % estimate model specified above
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
        % compute contrasts to get t-maps
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        
        switch iCon
            case 1
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'scenes > objects';
            case 2
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'scenes > scramble';
            otherwise
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'objects > scramble';
        end
    
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
        % go!
        spm_jobman('run_nogui', matlabbatch)
        clear matlabbatch;
    end
    
    %% 2nd-level GLM: Experimental Runs
    matlabbatch = [];
    numCons = 104;
    
    %{
        --> Why 104 contrasts? 
        1-100. per-image contrasts
        101. good > bad images
        102. faces > scenes
        103. natural > manmade
        104. foreground object > no foreground object
    %}
    
    for iCon = 1:numCons
        conFolderName = ['con' num2str(iCon)];
    
        % create folder for contrast if it doesn't exist
        if exist(fullfile(outPath, 'exp', conFolderName)) > 0
            rmdir(fullfile(outPath, 'exp', conFolderName), 's');
        end
    
        mkdir(fullfile(outPath, 'exp', conFolderName));
    
        % specify model
        matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(outPath, 'exp', conFolderName)};
    
        % get each subject's files for the contrast of interest
        conFiles = cell(length(subs), 1);
    
        for iSub = 1:subs
            % get subject ID
            if length(num2str(iSub)) < 2
                subID = ['sub-00' num2str(iSub) 'xxxx'];
            else
                subID = ['sub-0' num2str(iSub) 'xxxx'];
            end
    
            % ensure correct contrast ID (with correct amount of preceding 0's)
            if length(num2str(iCon)) == 1
                conID = ['000' num2str(iCon)];
            elseif length(num2str(iCon)) == 2
                conID = ['00' num2str(iCon)];
            else
                conID = ['0' num2str(iCon)];
            end
    
            % find contrast file
            conFiles{iSub,1} = fullfile( ...
                fmriPath, ...
                subID, ...
                sprintf( ...
                    'derivatives/exp_glm1_norm/con_%s.nii,1', conID ...
                ) ...
            );
        end
    
        % specify model
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = conFiles;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
        % estimate model specified above
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
        % compute contrasts to get t-maps
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
        switch iCon
            case num2cell(1:100)
                imgName = ['img' num2str(iCon)];
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = imgName;
            case 101
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'good > bad';
            case 102
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'faces > scenes';
            case 103
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'natural > manmade';
            otherwise
                matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'foreground > no foreground';
        end
    
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
        % go!
        spm_jobman('run_nogui', matlabbatch)
        clear matlabbatch;
    end
end

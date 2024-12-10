function level1_glm_norm_func(subs, runwise, nRuns)

%% make multiple condition files
sortRows = true;
includeTargets = false;
create_mcf_func(subs, sortRows, includeTargets)

%% Initialize SPM
matlabbatch = [];
spm('defaults', 'fmri');
spm_jobman('initcfg');

%% Define subjects and main path
mainPath = fullfile(pwd, '..');
fmriPath = fullfile(mainPath, 'sourcedata');

for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end

    % make `derivatives` sub-directory if it doesn't exist yet
    if ~exist(fullfile(mainPath, 'derivatives', subID), 'dir')
        mkdir(fullfile(mainPath, 'derivatives', subID));
    end

     %% Estimate model for localizer
    % make localizer sub-sub-directory if it doesn't exist yet
    if ~exist(fullfile(mainPath, 'derivatives', subID, 'loc_glm1_norm'), 'dir')
        mkdir(fullfile(mainPath, 'derivatives', subID, 'loc_glm1_norm'));
    end

    % model specification
    matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(mainPath, 'derivatives', subID, 'loc_glm1_norm')};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.85;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    % get each localizer scan's `.nii` file
    locFiles = dir(fullfile(mainPath, 'derivatives', subID, 'func', ...
        sprintf('wr%s%s_task-localizer_bold_*.nii', subID, 'xxxx')));

    locPaths = cell(1, length(locFiles));
    for i = 1:length(locFiles)
        locPaths{i} = [fullfile(locFiles(i).folder, locFiles(i).name), ',1'];
    end

    % get `mcf` file
    mcf = fullfile(mainPath, 'localizer', 'onsets',...
        sprintf('mcf_%s_localizer.mat', subID));

    % get motion regressors
    moRegs = fullfile(mainPath, 'derivatives', subID, 'func',...
        sprintf('rp_%s%s_task-localizer_bold_00001.txt', subID, 'xxxx'));

    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = locPaths';  % needs to be single-column cell array (hence transpose)
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {mcf};  % multiple-condition file goes here
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {moRegs};  % motion regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % contrasts
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'scenes > objects';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'scenes > scramble';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 0 -1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'objects > scramble';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    % go!
    spm_jobman('run', matlabbatch)
    clear matlabbatch;

    %% Estimate model for experimental runs

    % make experimental-run sub-sub-directory if it doesn't exist yet
    outputDir = fullfile(mainPath, 'derivatives', subID, 'exp_glm1_norm', filesep);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % for each run, get each scan's `.nii` file
    for run = 1:nRuns

        % make run directory if needed
        if runwise
            outputDir = fullfile(mainPath, 'derivatives', subID,...
                'exp_glm1_norm', sprintf('run-%02d', run), filesep);
            if ~exist(outputDir, 'dir')
                mkdir(outputDir);
            end
        end

        % do it before each run or only before the first one
        if run == 1 || runwise

            % clear specs
            matlabbatch = [];

            % model specification
            matlabbatch{1}.spm.stats.fmri_spec.dir = {outputDir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.85;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        end

        % get functional files
        expFiles = dir(fullfile(mainPath, 'derivatives', subID, 'func', ...
            sprintf('wr%s%s_task-scenes_run-%s_bold_*.nii',...
            subID, 'xxxx', num2str(run))));

        expPaths = cell(1, length(expFiles));
        for i = 1:length(expFiles)
            expPaths{i} = [fullfile(expFiles(i).folder, expFiles(i).name), ',1'];
        end
        
        % get `mcf` file
        mcf = fullfile(fmriPath, subID, 'beh', 'onsets', ...
            sprintf('mcf_%s_run-%s.mat', subID, num2str(run)));
        mcf_file = load(mcf);

        % get motion regressors
        moRegs = fullfile(mainPath, 'derivatives', subID, 'func',...
            sprintf('rp_%s%s_task-scenes_run-%s_bold_00001.txt',...
            subID, 'xxxx', num2str(run)));

        % specify session information
        if runwise
            sessNum = 1;
        else
            sessNum = run;
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).scans = expPaths';
        matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).multi = {mcf};
        matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).multi_reg = {moRegs};
        matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).hpf = 128;

        % do it after each run or only after the last one
        if run == nRuns || runwise

            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

            % model estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            % contrasts
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

            % per-image contrasts (only non-targets)
            totalNonTagertImg = 100;
            for img = 1:totalNonTagertImg
                matlabbatch{3}.spm.stats.con.consess{img}.tcon.name = mcf_file.names{img};
                contrastVec = [zeros(1, img-1), 1, zeros(1, totalNonTagertImg - img)];
                contrastVec = contrastVec - mean(contrastVec);
                matlabbatch{3}.spm.stats.con.consess{img}.tcon.weights = contrastVec;
                matlabbatch{3}.spm.stats.con.consess{img}.tcon.sessrep = 'repl';
            end

            % bathroom > kitchen images
            matlabbatch{3}.spm.stats.con.consess{img+1}.tcon.name = 'bathroom > kitchen';
            matlabbatch{3}.spm.stats.con.consess{img+1}.tcon.weights = [ones(1, 50) -(ones(1, 50))];
            matlabbatch{3}.spm.stats.con.consess{img+1}.tcon.sessrep = 'repl';

%             % target > non targett
%             matlabbatch{3}.spm.stats.con.consess{img+2}.tcon.name = 'target > nonTarget';
%             contrastVec = [zeros(1, totalNonTagertImg), ones(1, 16)];
%             contrastVec = contrastVec - mean(contrastVec);
%             matlabbatch{3}.spm.stats.con.consess{img+2}.tcon.weights = contrastVec;
%             matlabbatch{3}.spm.stats.con.consess{img+2}.tcon.sessrep = 'repl';

            % clear contrasts?
            matlabbatch{3}.spm.stats.con.delete = 0;  % no!

            % go!
            spm_jobman('run', matlabbatch)
            clear matlabbatch;
            disp(['Estimated (normalized) model for subject ', ...
                num2str(subs(iSub)), ' run ', num2str(run), '!']);

        end
    end
end
end

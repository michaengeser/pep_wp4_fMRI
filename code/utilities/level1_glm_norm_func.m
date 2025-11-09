function level1_glm_norm_func(cfg)

%
if ~isfield(cfg, 'hpf'); cfg.hpf = 128; end
if ~isfield(cfg, 'runwise'); cfg.runwise = true; end

%% make multiple condition files
cfg.sortRows = true;
create_mcf_func(cfg)

%% Initialize SPM
matlabbatch = [];
spm('defaults', 'fmri');
spm_jobman('initcfg');

%% Define subjects and main path
mainPath = fullfile(pwd, '..');
fmriPath = fullfile(mainPath, 'sourcedata');

for iSub = 1:length(cfg.subNums)

    % subject ID
    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

    % progress report
    disp(['Starting first level glm in SPM for subject ',...
        num2str(cfg.subNums(iSub))]);

    % make `derivatives` sub-directory if it doesn't exist yet
    if ~exist(fullfile(mainPath, 'derivatives', subID), 'dir')
        mkdir(fullfile(mainPath, 'derivatives', subID));
    end

    %% Estimate model for localizer
    % make a new localizer sub-sub-directory
    locOutputDir = fullfile(mainPath, 'derivatives', subID, 'loc_glm1_norm');
    skipLoc = false;
    if ~exist(locOutputDir, 'dir')
        mkdir(locOutputDir);
    else
        contrastFile = fullfile(locOutputDir, 'con_0004.nii');
        if exist(contrastFile, 'file')
            skipLoc = true;
            disp(['Data for subject ', num2str(cfg.subNums(iSub)), ' exists already'])
        else
            rmdir(locOutputDir, 's');
            pause(1); % Small delay before retry
            mkdir(locOutputDir);
        end
    end

    if ~skipLoc

        % model specification
        matlabbatch{1}.spm.stats.fmri_spec.dir = {locOutputDir};
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
        mcf_file = load(mcf);

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
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0, 1, -1, 0];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'scenes > scramble';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0, 1, 0, -1];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'objects > scramble';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0, 0, 1, -1];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'intact > scramble';
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1/3, 1/3, 1/3, -1];
        matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;

        % go!
        % repeat 3 times when crashing
        maxRetries = 3;
        for attempt = 1:maxRetries
            try
                spm_jobman('run', matlabbatch);
                break; % Exit loop if successful
            catch ME
                if attempt == maxRetries
                    warning('Failed on %s after %d attempts: %s', subID, maxRetries, ME.message);
                else
                    disp(['Retrying: Attempt ', num2str(attempt + 1)]);
                    cd(fullfile(pwd, '..'));
                    rmdir(locOutputDir, 's');
                    pause(1); % Small delay before retry
                    mkdir(locOutputDir);
                end
            end
        end

        clear matlabbatch;
    end

     %% Estimate model for experimental runs
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % define folder name
%     folderName = 'exp_glm1_norm';
% 
%     % make experimental-run sub-sub-directory if it doesn't exist yet
%     outputDir = fullfile(mainPath, 'derivatives', subID, folderName, filesep);
%     skipGLM = true;
%     if ~exist(outputDir, 'dir')
%         mkdir(outputDir);
%     else
%         contrastFile = fullfile(outputDir, 'con_0001.nii');
%         if exist(contrastFile, 'file')
%             skipGLM = true;
%             disp(['Data for subject ', num2str(cfg.subNums(iSub)), ' exists already'])
%         else
%             rmdir(outputDir, 's');
%             pause(1); % Small delay before retry
%             mkdir(outputDir);
%         end
%     end
% 
%     if ~skipGLM
%         % for each run, get each scan's `.nii` file
%         for run = 1:cfg.nRuns
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % make run directory if needed
%             if cfg.runwise
%                 outputDir = fullfile(mainPath, 'derivatives', subID,...
%                     folderName, sprintf('run-%02d', run), filesep);
%                 if ~exist(outputDir, 'dir')
%                     mkdir(outputDir);
%                 end
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             % do it before each run or only before the first one
%             if run == 1 || cfg.runwise
% 
%                 % clear specs
%                 matlabbatch = [];
% 
%                 % model specification
%                 matlabbatch{1}.spm.stats.fmri_spec.dir = {outputDir};
%                 matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
%                 matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.85;
%                 matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
%                 matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%             end
% 
%             % get functional files
%             expFiles = dir(fullfile(mainPath, 'derivatives', subID, 'func', ...
%                 sprintf('wr%s%s_task-scenes_run-%s_bold_*.nii',...
%                 subID, 'xxxx', num2str(run))));
% 
%             expPaths = cell(1, length(expFiles));
%             for i = 1:length(expFiles)
%                 expPaths{i} = [fullfile(expFiles(i).folder, expFiles(i).name), ',1'];
%             end
% 
%             % get `mcf` file
%             mcf = fullfile(fmriPath, subID, 'beh', 'onsets', ...
%                 sprintf('mcf_%s_run-%s.mat', subID, num2str(run)));
%             mcf_file = load(mcf);
% 
%             % get motion regressors
%             moRegs = fullfile(mainPath, 'derivatives', subID, 'func',...
%                 sprintf('rp_%s%s_task-scenes_run-%s_bold_00001.txt',...
%                 subID, 'xxxx', num2str(run)));
% 
%             % specify session information
%             if cfg.runwise
%                 sessNum = 1;
%             else
%                 sessNum = run;
%             end
%             matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).scans = expPaths';
%             matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
%             matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).multi = {mcf};
%             matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).regress = struct('name', {}, 'val', {});
%             matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).multi_reg = {moRegs};
%             matlabbatch{1}.spm.stats.fmri_spec.sess(sessNum).hpf = cfg.hpf;
% 
%             % do it after each run or only after the last one
%             if run == cfg.nRuns || cfg.runwise
% 
%                 matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
%                 matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
%                 matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
%                 matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
%                 matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
%                 matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
%                 matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
% 
%                 % model estimation
%                 matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%                 matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
%                 matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
% 
%                 % contrasts
%                 matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
% 
%                 % per-image contrasts
%                 runTrials = unique(cell2mat(mcf_file.trialIDs))';
%                 for itrialID = 1:length(runTrials)
%                     trialID = runTrials(itrialID);
%                     imageName = mcf_file.names(cell2mat(mcf_file.trialIDs) == trialID);
%                     matlabbatch{3}.spm.stats.con.consess{itrialID}.tcon.name = imageName{1};
%                     contrastVec = cell2mat(mcf_file.trialIDs) == trialID;
%                     contrastVec = contrastVec - mean(contrastVec);
%                     matlabbatch{3}.spm.stats.con.consess{itrialID}.tcon.weights = contrastVec;
%                     matlabbatch{3}.spm.stats.con.consess{itrialID}.tcon.sessrep = 'repl';
%                 end
% 
% 
%                 % clear contrasts?
%                 matlabbatch{3}.spm.stats.con.delete = 0;  % no!
% 
%                             % go!
%                             % repeat 3 times when crashing
%                             maxRetries = 3;
%                             for attempt = 1:maxRetries
%                                 try
%                                     spm_jobman('run', matlabbatch);
%                                     break; % Exit loop if successful
%                                 catch ME
%                                     if attempt == maxRetries
%                                         error('Failed after %d attempts: %s', maxRetries, ME.message);
%                                     else
%                                         disp(['Retrying: Attempt ', num2str(attempt)]);
%                                         pause(1); % Small delay before retry
%                                     end
%                                 end
%                             end
%                             disp(['Estimated (normalized) model for subject ', ...
%                                 num2str(cfg.subNums(iSub)), ' run ', num2str(run), '!']);
% 
% 
%                 % save trial IDs
%                 trialIDs = mcf_file.trialIDs;
%                 save(fullfile(outputDir, 'trialIDs.mat'), "trialIDs")
%             end
%         end
%     end
end
end

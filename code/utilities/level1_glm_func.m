function level1_glm_func(subs)    
    %% Initialize SPM
    matlabbatch = [];
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    
    %% Define subjects and main path
    subs = cellstr(string(subs));
    mainPath = '/Users/philippflieger/PhD/projects/fusion-study_scene-attractiveness/';
    fmriPath = fullfile(mainPath, 'sourcedata/fmri/bids');
    
    for sub = 1:length(subs)
        subID = ['sub-0', subs{sub}, 'xxxx'];
    
        % make `derivatives` sub-directory if it doesn't exist yet
        if ~exist(fullfile(fmriPath, subID, 'derivatives'))
            mkdir(fullfile(fmriPath, subID, 'derivatives'));
        end
    
        %% Estimate model for localizer
        % make localizer sub-sub-directory if it doesn't exist yet
        if ~exist(fullfile(fmriPath, subID, 'derivatives/loc_glm1'))
            mkdir(fullfile(fmriPath, subID, 'derivatives/loc_glm1'));
        end
    
        % model specification
        matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(fmriPath, subID, 'derivatives/loc_glm1')};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.85;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
        % get each localizer scan's `.nii` file
        locFiles = dir( ...
            fullfile( ...
                fmriPath, ...
                subID, ...
                'func', ...
                sprintf('r%s_task-localizer_bold_*.nii', subID) ...
            ) ...
        );
        
        locPaths = cell(1, length(locFiles));
        for i = 1:length(locFiles)
            locPaths{i} = [fullfile(locFiles(i).folder, locFiles(i).name), ',1'];
        end
    
        % get `mcf` file
        mcf = fullfile( ...
            mainPath, ...
            sprintf( ...
                'sourcedata/localizer/onsets/mcf_sub-%s_localizer.mat', ...
                subs{sub} ...
            ) ...
        );
    
        % get motion regressors
        moRegs = fullfile( ...
            fmriPath, ...
            sprintf( ...
                '%s/func/rp_%s_task-localizer_bold_00001.txt', ...
                subID, ...
                subID ...
            ) ...
        );
    
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
        spm_jobman('run_nogui', matlabbatch)
        clear matlabbatch;
    
        %% Estimate model for experimental runs
        matlabbatch = [];
    
        % make experimental-run sub-sub-directory if it doesn't exist yet
        if ~exist(fullfile(fmriPath, subID, 'derivatives/exp_glm1'))
            mkdir(fullfile(fmriPath, subID, 'derivatives/exp_glm1'));
        end
    
        % model specification
        matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(fmriPath, subID, 'derivatives/exp_glm1')};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.85;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
        % for each run, get each scan's `.nii` file
        for run = 1:7
            expFiles = dir( ...
                fullfile( ...
                    fmriPath, ...
                    subID, ...
                    'func', ...
                    sprintf('r%s_task-scenes_run-%s_bold_*.nii', subID, num2str(run)) ...
                ) ...
            );
            
            expPaths = cell(1, length(expFiles));
            for i = 1:length(expFiles)
                expPaths{i} = [fullfile(expFiles(i).folder, expFiles(i).name), ',1'];
            end
    
            % get `mcf` file
            subIDnoExt = sprintf('sub-%s', subs{sub});
            mcf = fullfile( ...
                mainPath, ...
                'sourcedata/beh', ...
                subIDnoExt, ...
                'onsets', ...
                sprintf('mcf_%s_run-%s.mat', subIDnoExt, num2str(run)) ...
            );
    
            % get motion regressors
            moRegs = fullfile( ...
                fmriPath, ...
                sprintf( ...
                    '%s/func/rp_%s_task-scenes_run-%s_bold_00001.txt', ...
                    subID, ...
                    subID, ...
                    num2str(run) ...
                ) ...
            );
    
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).scans = expPaths';
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi = {mcf};
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi_reg = {moRegs};
            matlabbatch{1}.spm.stats.fmri_spec.sess(run).hpf = 128;
        end
    
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
       
        % per-image contrasts
        for img = 1:100  
            matlabbatch{3}.spm.stats.con.consess{img}.tcon.name = sprintf('img%s', num2str(img));
            matlabbatch{3}.spm.stats.con.consess{img}.tcon.weights = [zeros(1, img-1) 1];
            matlabbatch{3}.spm.stats.con.consess{img}.tcon.sessrep = 'repl';
        end
    
        % good > bad images
        matlabbatch{3}.spm.stats.con.consess{img+1}.tcon.name = 'good > bad';
        matlabbatch{3}.spm.stats.con.consess{img+1}.tcon.weights = [ones(1, 50) -(ones(1, 50))];
        matlabbatch{3}.spm.stats.con.consess{img+1}.tcon.sessrep = 'repl';
        
        % faces > scenes
        matlabbatch{3}.spm.stats.con.consess{img+2}.tcon.name = 'faces > scenes';
        matlabbatch{3}.spm.stats.con.consess{img+2}.tcon.weights = convertWeights('faceWeights');
        matlabbatch{3}.spm.stats.con.consess{img+2}.tcon.sessrep = 'repl';
    
        % natural > manmade
        matlabbatch{3}.spm.stats.con.consess{img+3}.tcon.name = 'natural > manmade';
        matlabbatch{3}.spm.stats.con.consess{img+3}.tcon.weights = convertWeights('naturalWeights');
        matlabbatch{3}.spm.stats.con.consess{img+3}.tcon.sessrep = 'repl';
    
        % foreground object > no foreground object
        matlabbatch{3}.spm.stats.con.consess{img+4}.tcon.name = 'foreground > no foreground';
        matlabbatch{3}.spm.stats.con.consess{img+4}.tcon.weights = convertWeights('foregroundWeights');
        matlabbatch{3}.spm.stats.con.consess{img+4}.tcon.sessrep = 'repl';
    
        % clear contrasts?
        matlabbatch{3}.spm.stats.con.delete = 0;  % no!
        
        % go!
        spm_jobman('run_nogui', matlabbatch)
        clear matlabbatch;
        disp(['Estimated model for subject ', char(subs(sub)), '!']);
    end
end
        
function updatedWeights = convertWeights(weightsFile)
    arguments
        weightsFile {mustBeText}
    end

    w = load(sprintf('%s.mat', weightsFile));
    w = w.stim_category';
    w(w==2) = -1;
    updatedWeights = w;
end

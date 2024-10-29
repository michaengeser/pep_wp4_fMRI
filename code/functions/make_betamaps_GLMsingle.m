%{
    This script is used to estimate beta maps for each trial of an fMRI
    experiment. Originally written by Lixiang Chen, adapted by Philipp
    Flieger (September 2024).
%}

clear; close all; clc;

% add toolboxes
addpath(genpath('GLMsingle/matlab/'));
addpath(genpath('GLMsingle/matlab/utilities/'));
addpath(genpath('fracridge/'));
addpath(genpath('NIfTI_20140122/'));

% set input and output paths
mainPath = '/Users/philippflieger/PhD/projects/fusion-study_scene-attractiveness/';
fmriPath = fullfile(mainPath, 'sourcedata/fmri/bids');
behPath = fullfile(mainPath, 'sourcedata/beh');

% set parameters
subs = [1:4 6:30]; % subject numbers
nRuns = 7; % number of runs
nConds = 1; % number of conditions
nImgs = 100; % number of scenes
TR = 1.85; % TR = 1.85s
stimDur = 1.45; % duration of stimulus
filePrefix = 'wr';  % file prefix for realigned, normalized volumes

%% Estimate betas
for iSub = 1:length(subs)
    subID = sprintf('sub-%03dxxxx', subs(iSub));

    if ~exist( ...
        fullfile(fmriPath, subID, 'derivatives/GLMsingle_betas/beta.nii'), ...
        'file' ...
    )
        fprintf('Creating beta map for subject %d...\n', subs(iSub));

	    % prepare data
        data = {}; % fMRI data
        design = {}; % event information
        motion = {}; % motion information
    
        for iRun = 1:nRuns
            runID = sprintf('run-%s', num2str(iRun));
    
            % create 4D file from realigned, normalized 3D volumes (if needed)
            realignedFile = [
                filePrefix subID '_task-scenes_run-' num2str(iRun) '_bold.nii' ...
            ];
    
            if ~exist(fullfile(fmriPath, subID, 'func', realignedFile), 'file')
                % get 3D files for run in question
                files3D = dir( ...
                    fullfile( ...
                        fmriPath, subID, 'func', ...
                        sprintf( ...
                            'wr%s_task-scenes_run-%s_bold_*.nii', ...
                            subID, num2str(iRun) ...
                        ) ...
                    ) ...
                );
                paths3D = cell(1, length(files3D));
    
                for iPath = 1:length(files3D)
                    paths3D{iPath} = [ ...
                        fullfile(files3D(iPath).folder, files3D(iPath).name) ',1' ...
                    ];
                end
    
                % init. SPM
                spm('defaults', 'fmri');
                spm_jobman('initcfg');
                matlabbatch = [];
    
                % define 3D-to-4D file conversion
                matlabbatch{1}.spm.util.cat.vols = paths3D';
                matlabbatch{1}.spm.util.cat.name = realignedFile;
                matlabbatch{1}.spm.util.cat.dtype = 4;
                matlabbatch{1}.spm.util.cat.RT = TR;
    
                % go!
                spm_jobman('run_nogui', matlabbatch);
                clear matlabbatch;
            end
    
		    % load pre-processed 4D fMRI data
            funcData = fullfile(fmriPath, subID, 'func', realignedFile);
		    runData = load_untouch_nii(funcData);
            data{iRun} = runData.img;
    
            % load `mcf` file (already contains onsets for each stimulus)
            mcf = load( ...
                fullfile( ...
                    behPath, sprintf('sub-%02d', subs(iSub)), 'onsets', ...
                    sprintf('mcf_sub-%02d_run-%d.mat', subs(iSub), iRun) ...
                ) ...
            );
    
            % create design matrix indexing the volume each stimulus occurs in
            nVols3D = size(runData.img, 4);
            designMatrix = zeros(nVols3D, nConds * nImgs);  % initialize
            onsetVolsByImg = ceil(cell2mat(mcf.onsets) / TR)';
    
            % fill design matrix
            for iOnset = 1:length(onsetVolsByImg)
                designMatrix(onsetVolsByImg(iOnset), iOnset) = 1;
            end
    
            % update `data` and `design` cell arrays
            data{iRun} = single(runData.img);
            design{iRun} = sparse(designMatrix);
		    
		    % load motion regressors
		    moRegs = fullfile( ...
                fmriPath, subID, 'func', ...
                ['rp_' subID '_task-scenes_' runID '_bold_00001.txt'] ...
            );
            motion{iRun} = load(moRegs);
        end
        
        %% Estimate beta using GLMsingle
        % create output folder
        betaOutDir = fullfile(fmriPath, subID, 'derivatives/GLMsingle_betas');
        if ~exist(betaOutDir, 'dir'); mkdir(betaOutDir); end
    
        % estimate beta and save it in output directory
        opt.wantmemoryoutputs = [1, 1, 1, 1]; 
        opt.extraregressors = motion;
        GLMestimatesingletrial(design, data, stimDur, TR, betaOutDir, opt);
        
        %% convert beta mat to nii
        betaModel = 'TYPED_FITHRF_GLMDENOISE_RR'; % FITHRF_GLMdenoise_RR model
            
        % load beta model
        load(fullfile(betaOutDir, betaModel), 'modelmd');
        beta = modelmd;
        
        % load event information  
        event = [];
    
        for iRun = 1:nRuns
            behData = dir( ...
                fullfile( ...
                    behPath, sprintf('sub-%02d', subs(iSub)), ...
                    sprintf('subj%d_run%d*.mat', subs(iSub), iRun) ...
                ) ...
            );
            load(fullfile(behData.folder, behData.name), 'dat');
            
            % get order of stimulus presentation
            cond = 1;
    
            % ensure correct image indices
            scene = zeros(nImgs, 1);
    
            for iRow = 1:length(dat.randomTrials)
                if dat.randomTrials(iRow,1) == 1
                    scene(iRow) = dat.randomTrials(iRow, 2);
                else
                    scene(iRow) = dat.randomTrials(iRow, 2) + 50;
                end
            end
            
            % add everything to `event` array
            singleRunEventArray = [];
    
            for iScene = 1:length(scene)
                singleRunEventArray = [ ...
                    singleRunEventArray; [iRun cond scene(iScene)] ...
                ];
            end
    
            event = [event; singleRunEventArray];
        end
    
        % sort event
        event = [event (1:size(event, 1))'];
        event = sortrows(event, [1 2 3]);  % run * cond * scene
        
        % sort beta
        beta = beta(:, :, :, event(:,4));  % run * cond * scene
        
        % save beta as nifti file
        ref = fullfile(fmriPath, subID, 'func', realignedFile);
        data = load_untouch_nii(ref);
        data.img = beta;
        data.hdr.dime.datatype = 16;
        data.hdr.dime.dim(5) = size(beta, 4);
        save_untouch_nii(data, fullfile(betaOutDir, 'beta.nii'));
    else
        fprintf('Beta map for subject %d already exists.\n\n', subs(iSub));
    end
end

%{
    This script is used to estimate beta maps for each trial of an fMRI
    experiment. Originally written by Lixiang Chen, adapted by Philipp
    Flieger (September 2024).
%}

function make_betamaps_GLMsingle(subs, nRuns, nTrials)

% add toolboxes
addpath(genpath(fullfile(pwd,'utilities','GLMsingle','matlab')));
addpath(genpath(fullfile(pwd,'utilities','GLMsingle','matlab','utilities',filesep)));
addpath(genpath(fullfile(pwd,'utilities','GLMsingle','matlab','fracridge',filesep)));
addpath(genpath(fullfile(pwd,'utilities','NIfTI_20140122',filesep)));

% set input and output paths
mainPath = fullfile(pwd,'..');
dataPath = fullfile(mainPath, 'sourcedata');
dervPath = fullfile(mainPath, 'derivatives');

% set parameters
nImgs = nTrials; % number of total scenes
TR = 1.85; % TR = 1.85s
stimDur = .25; % duration of stimulus
filePrefix = 'wr';  % file prefix for realigned, normalized volumes

%% make multiple condition files
sortRows = true;
includeTargets = true;
create_mcf_func(subs, sortRows, includeTargets)

%% Estimate betas
for iSub = 1:length(subs)

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
        subIDx = ['sub-00', num2str(subs(iSub)), 'xxxx'];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
        subIDx = ['sub-0', num2str(subs(iSub)), 'xxxx'];
    end

    if ~exist(fullfile(dervPath, subID, 'GLMsingle_betas', 'beta.nii'), 'file')

        disp(char(datetime));
        fprintf('Creating beta map for subject %d...\n', subs(iSub));

        % prepare data
        data = {}; % fMRI data
        design = {}; % event information
        motion = {}; % motion information

        for iRun = 1:nRuns
           fprintf('Processing run %d...\n', iRun);

            % create 4D file from realigned, normalized 3D volumes (if needed)
            realignedFile = [filePrefix, subIDx, '_task-scenes_run-',...
                num2str(iRun), '_bold.nii'];

            if ~exist(fullfile(dervPath, subID, 'func', realignedFile), 'file')
                % get 3D files for run in question
                files3D = dir(fullfile(dervPath, subID, 'func', ...
                    sprintf('wr%s_task-scenes_run-%s_bold_*.nii', ...
                    subIDx, num2str(iRun))));

                paths3D = cell(1, length(files3D));

                for iPath = 1:length(files3D)
                    paths3D{iPath} = [fullfile(files3D(iPath).folder,...
                        files3D(iPath).name) ',1'];
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
            funcData = fullfile(dervPath, subID, 'func', realignedFile);
            runData = load_untouch_nii(funcData);
            data{iRun} = runData.img;

            % load `mcf` file (already contains onsets for each stimulus)
            mcf = load(fullfile(dataPath, subID, 'beh','onsets', ...
                sprintf('mcf_%s_run-%d.mat', subID, iRun)));

            % create design matrix indexing the volume each stimulus occurs in
            nVols3D = size(runData.img, 4);
            designMatrix = zeros(nVols3D, nImgs);  % initialize
            onsetVolsByImg = ceil(cell2mat(mcf.onsets) / TR)';

            % fill design matrix
            nRegressors = nImgs; % includes only non-targets
            % nRegressors = length(onsetVolsByImg); % includes also targets
            for iOnset = 1:nRegressors
                designMatrix(onsetVolsByImg(iOnset), iOnset) = 1;
            end

            % update `data` and `design` cell arrays
            data{iRun} = single(runData.img);
            design{iRun} = sparse(designMatrix);

            % load motion regressors
            moRegs = fullfile(dervPath, subID, 'func', ...
                ['rp_', subIDx, '_task-scenes_run-', num2str(iRun), '_bold_00001.txt']);
            motion{iRun} = load(moRegs);
        end

        %% Estimate beta using GLMsingle
        % create output folder
        betaOutDir = fullfile(dervPath, subID, 'GLMsingle_betas');
        if ~exist(betaOutDir, 'dir'); mkdir(betaOutDir); end

        % estimate beta and save it in output directory
        opt.wantmemoryoutputs = [0, 0, 0, 1];
        opt.extraregressors = motion;
        GLMestimatesingletrial(design, data, stimDur, TR, betaOutDir, opt);

        %% convert beta mat to nii
        betaModel = 'TYPED_FITHRF_GLMDENOISE_RR'; % FITHRF_GLMdenoise_RR model

        % load beta model
        load(fullfile(betaOutDir, betaModel), 'modelmd');
        beta = modelmd;

        % sorting beta
        % load event information
        event = [];

        for iRun = 1:nRuns
            if iRun < 10
                runStr = ['0', num2str(iRun)];
            else
                runStr = num2str(iRun);
            end

            % get log file 
            dat = readtable(fullfile(dataPath, subID, 'beh', ...
                sprintf('%s_task-main_run-%s_events.tsv', subID, runStr)),...
                'FileType', 'text', 'Delimiter', '\t');

            % remove target trials 
            dat = dat(~strcmp(dat.category, 'livingroom'), :);

            % ensure correct image indices
            scene = dat.texture' - 10;
            scene = scene(scene <= nTrials);

            % get condition
            cond = zeros(1, nTrials);
            for tr = 1:nTrials
                if strcmp(dat.category{tr}, 'bathroom')
                    cond(tr) = 1;
                else
                    cond(tr) = 2;
                end
            end

            % add everything to `event` array
            singleRunEventArray = [];

            for iScene = 1:length(scene)
                singleRunEventArray = [ ...
                    singleRunEventArray; [iRun, cond(iScene), scene(iScene)]];
            end

            event = [event; singleRunEventArray];
        end

        % sort event
        event = [event (1:size(event, 1))'];
        event = sortrows(event, [1 2 3]);  % run * cond * scene

        % sort beta
        beta = beta(:, :, :, event(:,4));  % run * cond * scene

        % save beta as nifti file
        ref = fullfile(dervPath, subID, 'func', realignedFile);
        data = load_untouch_nii(ref);
        data.img = beta;
        data.hdr.dime.datatype = 16;
        data.hdr.dime.dim(5) = size(beta, 4);
        save_untouch_nii(data, fullfile(betaOutDir, 'beta.nii'));
    else
        fprintf('Beta map for subject %d already exists.\n\n', subs(iSub));
    end
end
end
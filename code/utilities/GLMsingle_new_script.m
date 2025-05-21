function GLMsingle_new_script(cfg)

nRuns = 12;
tr = 1.85;
nTrials = 100;
nVols = 188;
stimdur = 0.25;
showReliability = false;

%% make multiple condition files
cfg.sortRows = true;
create_mcf_func(cfg)

%% Define subjects and main path
mainPath = fullfile(pwd, '..');
fmriPath = fullfile(mainPath, 'sourcedata');

%% get data

for iSub = 1:length(cfg.subNums)

    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

    % make `derivatives` sub-directory if it doesn't exist yet
    if ~exist(fullfile(mainPath, 'derivatives', subID), 'dir')
        mkdir(fullfile(mainPath, 'derivatives', subID));
    end

    % check if poutput exists already 
    outputdir = fullfile(mainPath, 'derivatives', subID, 'GLMsingleEstimates');
    if exist(outputdir, 'dir')
        outputFile = fullfile(outputdir, 'GLMsingle_betas.nii');
        if exist(outputFile, 'file')
            disp(['Data for subject ', num2str(cfg.subNums(iSub)), ' exists already'])
            continue
        end 
    else
        mkdir(outputdir);
    end

    %% Loop through runs

    % for each run, get each scan's `.nii` file
    data = cell(1, nRuns);
    design = cell(1, nRuns);
    trialIDs = [];
    %motion = cell(1, nRuns); % not  necessary

    for iRun = 1:nRuns

        % get functional files
        funcFile = fullfile(mainPath, 'derivatives', subID, 'func', ...
            sprintf('wr%s%s_task-scenes_run-%s_bold.nii',...
            subID, 'xxxx', num2str(iRun)));

        %% check if 4D file exists
        if ~exist(funcFile, 'file')

            % get 3D files
            funcFiles = dir(strrep(funcFile, 'bold', 'bold_0*'));

            if height(funcFiles) ~= nVols
                warning(['Number of volumnes is not ', num2str(nVols)])
            end

            % convert 3D files to 4D file
            disp(['Convert 3D files to 4D files for run ', num2str(iRun)])

            fileCell = cell(height(funcFiles), 1);
            for iVol = 1:height(funcFiles)

                fileCell{iVol, 1} = fullfile(funcFiles(iVol).folder,...
                    [funcFiles(iVol).name, ',1']);

            end
            % init. SPM
            spm('defaults', 'fmri');
            spm_jobman('initcfg');
            matlabbatch = [];
            
            % run spm
            matlabbatch{1}.spm.util.cat.vols = fileCell;
            matlabbatch{1}.spm.util.cat.name = funcFile;
            matlabbatch{1}.spm.util.cat.dtype = 4;
            matlabbatch{1}.spm.util.cat.RT = tr;
            spm_jobman('run',matlabbatch)
        end

        %% load data
        v = load_untouch_nii(funcFile);
        data{iRun} = single(v.img);

        % get `mcf` file
        mcf = fullfile(fmriPath, subID, 'beh', 'onsets', ...
            sprintf('mcf_%s_run-%s.mat', subID, num2str(iRun)));
        mcf_file = load(mcf);

        % add trial IDs
        trialIDs = [trialIDs; cell2mat(mcf_file.trialIDs),...
            ones(size(mcf_file.trialIDs))*iRun,...
            cell2mat(mcf_file.onsets)];

        % check if number of conditions match onsets
        if ~numel(mcf_file.onsets) == nTrials
            warning('Onsets does not match number of conditions')
        end

        %% make design matrix
        design{iRun} = zeros(nVols, nTrials);

        % set onset in correct TR
        for iTrials=1:size(design{iRun}, 2)
            if ismember(iTrials, cell2mat(mcf_file.trialIDs))

                % get idx of condition column
                condIdx = find(iTrials == cell2mat(mcf_file.trialIDs));

                % set correct tr to 1
                design{iRun}(round(cell2mat(mcf_file.onsets(condIdx)) / tr) + 1,...
                    iTrials) = 1;
            end
        end

        %     % get motion regressors
        %     moRegs = fullfile(mainPath, 'derivatives', subID, 'func',...
        %         sprintf('rp_%s%s_task-scenes_run-%s_bold_00001.txt',...
        %         subID, 'xxxx', num2str(iRun)));
        %     motion{iRun} = load(moRegs);

    end

    %% run GLM single

    % set options
    opt = struct('wantmemoryoutputs',[0 0 0 1]); % only type d model should be written to memory
    opt.wantfileoutputs = [0 0 0 1];
    %opt.extraregressors = motion;

    % get model estimate
    [results, ~] = GLMestimatesingletrial(design, data, stimdur, tr, outputdir, opt);
    model = results{4}.modelmd;


    %% save results 

    % sort and save trials IDs
    trialIDs = sortrows(trialIDs, [2, 3]);
    save(fullfile(outputdir, 'trialIDs.mat'), 'trialIDs')

    %     % save model estimates
    %     V.fname = 'GLMsingle_betas.nii';
    %     V.dim = size(model);
    %     V.dt = [16, 0]; % Data type (e.g., 16 = float32)
    %     V.mat = eye(4); % Identity matrix for affine transformation
    %     V.descrip = 'GLMsingle TypeD model estimates NIfTI file';

    % Write NIfTI file
    v.img = model;
    v.hdr.dime.datatype = 16;
    v.hdr.dime.dim(2:5) = size(model);
    save_untouch_nii(v, fullfile(outputdir, 'GLMsingle_betas.nii'));


    %% get reliability

    if showReliability
        % Consolidate design matrices
        designALL = cat(1,design{:});

        % Construct a vector containing 1-indexed condition numbers in
        % chronological order.

        corder = [];
        for p=1:size(designALL,1)
            if any(designALL(p,:))
                corder = [corder find(designALL(p,:))];
            end
        end

        % Get indices of repeated stimuli (not targets, columns are conditions, 
        % rows are runs/repetitions)
        repindices = [];
        for p=1:nTrials
            temp = find(corder==p);
            repindices = cat(2,repindices,temp');  % note that for images with 3 presentations, we are simply ignoring the third trial
        end

        % get voxel reliabilities
        pairs = nchoosek(1:10,2);
        vox_reliabilities = cell(1,length(pairs));

        % For each pair of trials
        for p = 1 : length(pairs)

            % Get the GLM betas
            betas = model(:,:,:,repindices);  % use indexing to pull out the trials we want
            betas_reshaped = reshape(betas, size(betas,1), size(betas,2), size(betas, 3),...
                nRuns, []);  % reshape to X x Y x Z x nRuns x nTrials

            % compute reliabilities using an efficient (vectorized) utility
            % function
            vox_reliabilities{p} = calccorrelation(...
                betas_reshaped(:, :, :, pairs(p, 1), :), ...
                betas_reshaped(:, :, :, pairs(p, 2), :), 5);

            % Note that calccorrelation.m is a utility function that computes
            % correlations in a vectorized fashion (for optimal speed).

        end

        % average across all pairs
        sumMatrix = zeros(size(vox_reliabilities{1}));
        for i = 1:numel(vox_reliabilities)
            sumMatrix = sumMatrix + vox_reliabilities{i};
        end
        averageMatrix = sumMatrix / numel(vox_reliabilities);


        % get visual cortex ROI
        roiDir = fullfile(pwd, '..', 'MNI_ROIs', 'wvisualCortex.nii');
        vROI = load_untouch_nii(roiDir);
        ROI = double(vROI.img);

        % display reliability
        disp(['The median reliability of voxels in the viusal cortex across repetitions is: ', ...
            num2str(nanmedian(averageMatrix(ROI==1)))])
    end
end

end
function GLMsingle_new_script(subs)

nRuns = 10;
tr = 1.85;
nTrials = 100;
nVols = 152;
stimdur = 0.25;
showReliability = true;
nTargets = 50;


%% make multiple condition files
sortRows = true;
includeTargets = true;
create_mcf_func(subs, sortRows, includeTargets)


%% Define subjects and main path
mainPath = fullfile(pwd, '..');
fmriPath = fullfile(mainPath, 'sourcedata');



%% get data

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

    % make make output fodler doesn't exist yet
    outputdir = fullfile(mainPath, 'derivatives', subID, 'GLMsingleEstimates');
    if ~exist(outputdir, 'dir')
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

        if ~exist(funcFile, 'file')
            error('make 4D file')
        end

        % load data
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

        % make design matrix
        if includeTargets
            design{iRun} = zeros(nVols, nTrials + nTargets);
        else
            design{iRun} = zeros(nVols, nTrials);
        end

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
        v = load_untouch_nii(roiDir);
        ROI = double(v.img);

        % display reliability
        disp(['The median reliability of voxels in the viusal cortex across repetitions is: ', ...
            num2str(nanmedian(averageMatrix(ROI==1)))])
    end
end

end
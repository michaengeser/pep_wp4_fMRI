function averageVoxelsInROI(cfg)

% evaluate input
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 12; end
if ~isfield(cfg, 'skipIfExists'); cfg.skipIfExists = false; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85;end
if ~isfield(cfg, 'nVols'); cfg.nVols = 188;end
if ~isfield(cfg, 'nTrials'); cfg.nTrials = 100;end
if ~isfield(cfg, 'nTargets'); cfg.nTargets = 15; end
if ~isfield(cfg, 'halfTrials'); cfg.halfTrials =  cfg.nTrials / 2;end
if ~isfield(cfg, 'TRstartBuffer'); cfg.TRstartBuffer = 3;end
if ~isfield(cfg, 'TRendBuffer'); cfg.TRendBuffer = 3;end
if ~isfield(cfg, 'targetDelay'); cfg.targetDelay = 3;end % in secs
if ~isfield(cfg, 'cutTargets'); cfg.cutTargets = true;end % in secs
if ~isfield(cfg, 'smoothing'); cfg.smoothing = false;end
if ~isfield(cfg, 'smoothKernel'); cfg.smoothKernel = 12; end
if ~isfield(cfg, 'saveWholeBrain'); cfg.saveWholeBrain = false;end
if ~isfield(cfg, 'rois'); cfg.rois = {'wV1.nii', 'wV2.nii',...
        'wLOC.nii', 'wPPA.nii',...
        'wTOS.nii', 'wRSC.nii'};
end

% init variables
subs = cfg.subNums;

%% make multiple condition files
cfg.sortRows = true;
cfg.includeTargets = true;
create_mcf_func(cfg)

%% load ROI mask
nmasks = numel(cfg.rois);

%% get data

for iSub = 1:length(subs)

    % runtime control
    disp(['Subject: ', num2str(subs(iSub))])

    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
    subID2 = strrep(subID, '-', '');

    % make `derivatives` sub-directory if it doesn't exist yet
    if ~exist(fullfile(cfg.outputPath, subID), 'dir')
        mkdir(fullfile(cfg.outputPath, subID));
    end

    % make make output folder doesn't exist yet
    if cfg.cutTargets
        outputdir = fullfile(cfg.outputPath, subID, 'timecourses');
    else
        outputdir = fullfile(cfg.outputPath, subID, 'timecourses_with_targets');
    end
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end

    %% Loop through runs

    % for each run, get each scan's `.nii` file
    data = cell(1, cfg.nRuns);

    for iRun = 1:cfg.nRuns

        % get category of the run
        if mod(iRun, 2) == 1
            currentCat = 'bathroom';
        else
            currentCat = 'kitchen';
        end

        % get functional files
        if cfg.smoothing
            if cfg.smoothKernel == 12
                funcFile = fullfile(cfg.outputPath, subID, 'func', ...
                    sprintf('s12wr%s%s_task-scenes_run-%s_bold.nii',...
                    subID, 'xxxx', num2str(iRun)));
            else
                funcFile = fullfile(cfg.outputPath, subID, 'func', ...
                    sprintf('swr%s%s_task-scenes_run-%s_bold.nii',...
                    subID, 'xxxx', num2str(iRun)));
            end
        else
            funcFile = fullfile(cfg.outputPath, subID, 'func', ...
                sprintf('wr%s%s_task-scenes_run-%s_bold.nii',...
                subID, 'xxxx', num2str(iRun)));
        end

        %% check if 4D file exists
        if ~exist(funcFile, 'file')

            % get 3D files
            funcFiles = dir(strrep(funcFile, 'bold', 'bold_0*'));

            if height(funcFiles) ~= cfg.nVols
                warning(['Number of volumnes is not ', num2str(cfg.nVols)])
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
            matlabbatch{1}.spm.util.cat.RT = cfg.tr;

            % repeat 3 times when crashing
            maxRetries = 3;
            for attempt = 1:maxRetries
                try
                    spm_jobman('run', matlabbatch);
                    break; % Exit loop if successful
                catch ME
                    if attempt == maxRetries
                        error('Failed after %d attempts: %s', maxRetries, ME.message);
                    else
                        disp(['Retrying: Attempt ', num2str(attempt)]);
                        pause(1); % Small delay before retry
                    end
                end
            end
        end

        %% load data and trial information

        % get `mcf` file from first subject (identical timing for all
        % subjects)
        mcf = fullfile(cfg.sourcedataPath, 'sub-101', 'beh', 'onsets', ...
            sprintf('mcf_sub-101_run-%s.mat', num2str(iRun)));
        mcf_file = load(mcf);

        % get trial IDs
        trialIDs = [cell2mat(mcf_file.trialIDs),...
            ones(size(mcf_file.trialIDs))*iRun,...
            cell2mat(mcf_file.onsets)];

        % sort trial IDs and add TR
        trialIDs = sortrows(trialIDs, 3);
        trialIDs(:, 4) = round(trialIDs(:, 3) / cfg.tr) + 1;

        % get target positions
        load(fullfile(pwd, 'utilities', 'targets.mat'), 'targetStruct')
        targetTrials = targetStruct(iRun).trialNum;

        % define vector which values to cut out
        includedTRs = ones(1, cfg.nVols);
        if cfg.cutTargets
            for iTr = 1:length(targetTrials)
                targetImg = trialIDs(targetTrials(iTr), 3) + cfg.targetDelay + cfg.stimDur + cfg.iti; % cut from the end of the trial
                targetImg = round(targetImg / cfg.tr) + 1; % convert to TR
                if targetTrials(iTr)  == cfg.nTrials % if target was last image
                    nextImg = cfg.nVols; % cut until the end
                else
                    nextImg = trialIDs(targetTrials(iTr) + 1, 3) + cfg.targetDelay;
                    nextImg = round(nextImg / cfg.tr) + 1; % convert to TR
                end
                includedTRs(targetImg:nextImg) = 0;
            end
        end

        % cut first and last TRs
        includedTRs(1:cfg.TRstartBuffer) = 0;
        includedTRs(end-cfg.TRendBuffer+1:end) = 0;

        % get data
        v = load_untouch_nii(funcFile);
        runData = single(v.img);
        data{iRun} = runData(:, :, :, logical(includedTRs));

        %% average across voxels for each mask
        for j=1:nmasks

            % get mask name
            mask_label=cfg.rois{j};
            mask_label_short = split(mask_label, '.');
            mask_label_short = mask_label_short{1};
            mask_label_short = mask_label_short(2:end);

            if cfg.smoothing
                if cfg.smoothKernel == 12
                    fileName = ['smoothed12_mean_timecourse_', currentCat, ...
                        '_', mask_label_short, ...
                        '_run_', num2str(iRun), '.mat'];
                else
                    fileName = ['smoothed_mean_timecourse_', currentCat, ...
                        '_', mask_label_short, ...
                        '_run_', num2str(iRun), '.mat'];

                end

               
            else
                fileName = ['mean_timecourse_', currentCat, ...
                    '_', mask_label_short, ...
                    '_run_', num2str(iRun), '.mat'];
            end

            % check if file exists already and skip if the case
            if cfg.skipIfExists

                if exist(fullfile(outputdir, fileName), 'file')
                    disp(['file for ', mask_label_short, ' and ', num2str(iRun), ...
                        ' exists already'])
                    continue
                end
            end

            % check if functional or anatomical ROI
            if ismember(mask_label_short, {'PPA', 'TOS', 'RSC', 'LOC', 'LPFC'})
                mask_fn=fullfile(pwd, '..', 'MNI_ROIs', 'func_ROIs', subID,...
                    [mask_label_short, '_funcROI.nii']);
            else
                mask_fn=fullfile(pwd, '..', 'MNI_ROIs', [char(mask_label)]);
            end

            % get mask
            mask = load_untouch_nii(mask_fn);
            newMaskImg =double(mask.img);
            if max(max(max(double(newMaskImg)))) > 1
                newMaskImg = newMaskImg/max(max(max(newMaskImg)));
            end

            % make 4D mask
            currentROI = newMaskImg;
            currentROI4D = repmat(currentROI, ...
                [1, 1, 1, size(data{iRun}, 4)]);

            % filter for data and take mean
            currentImage = data{iRun};
            filtered = currentImage(currentROI4D == 1);
            ROIImage = reshape(filtered, [], size(data{iRun}, 4));
            meanData = mean(ROIImage, 1);

            % save data
            saveData = meanData;
            save(fullfile(outputdir, fileName), 'saveData')
        end % rois

        if cfg.saveWholeBrain && cfg.smoothing

            % save whole brain data
            fileName = ['whole_brain_timecourse_', currentCat, ...
                '_run_', num2str(iRun), '.mat'];

            % check if file exists already and skip if the case
            if cfg.skipIfExists

                if exist(fullfile(outputdir, fileName), 'file')
                    disp(['file for whole brain and ', num2str(iRun), ...
                        ' exists already'])
                    continue
                end
            end

            % saving
            saveData = data{iRun};
            save(fullfile(outputdir, fileName), 'saveData')
        end
    end % runs
end % subjects
end
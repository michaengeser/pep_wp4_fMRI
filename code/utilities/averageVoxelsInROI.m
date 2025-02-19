function averageVoxelsInROI(cfg)

% evaluate input
if ~isfield(cfg, 'nRuns'); cfg.nRuns = 10; end
if ~isfield(cfg, 'skipIfExists'); cfg.skipIfExists = true; end
if ~isfield(cfg, 'tr'); cfg.tr = 1.85;end
if ~isfield(cfg, 'nVols'); cfg.nVols = 152;end
if ~isfield(cfg, 'nTrials'); cfg.nTrials = 100;end
if ~isfield(cfg, 'nTargets'); cfg.nTargets = 16; end
if ~isfield(cfg, 'halfTrials'); cfg.halfTrials =  cfg.nTrials / 2;end
if ~isfield(cfg, 'TRstartBuffer'); cfg.TRstartBuffer = 0;end
if ~isfield(cfg, 'TRendBuffer'); cfg.TRendBuffer = 4;end
if ~isfield(cfg, 'smoothing'); cfg.smoothing = false;end
if ~isfield(cfg, 'saveWholeBrain'); cfg.saveWholeBrain = false;end
if ~isfield(cfg, 'rois'); cfg.rois = {'wV1.nii', 'wV2.nii',...
        'wLOC.nii', 'wPPA.nii',...
        'wTOS.nii', 'wRSC.nii'};
end

% init variables
subs = cfg.subNums;

%% make multiple condition files
sortRows = true;
includeTargets = true;
create_mcf_func(subs, sortRows, includeTargets)

%% Define subjects and main path
mainPath = fullfile(pwd, '..');
fmriPath = fullfile(mainPath, 'sourcedata');

%% load ROI mask
nmasks = numel(cfg.rois);

%% get data

for iSub = 1:length(subs)


    % runtime control
    disp(['Subject: ', num2str(subs(iSub))])

    if subs(iSub) < 10
        subID = ['sub-00', num2str(subs(iSub))];
    elseif subs(iSub) < 100
        subID = ['sub-0', num2str(subs(iSub))];
    end
    subID2 = strrep(subID, '-', '');

    % make `derivatives` sub-directory if it doesn't exist yet
    if ~exist(fullfile(mainPath, 'derivatives', subID), 'dir')
        mkdir(fullfile(mainPath, 'derivatives', subID));
    end

    % make make output folder doesn't exist yet
    outputdir = fullfile(mainPath, 'derivatives', subID, 'timecourses');
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end

    %% Loop through runs

    % for each run, get each scan's `.nii` file
    data = cell(1, cfg.nRuns);
    bathroomData = data;
    kitchenData = data;

    for iRun = 1:cfg.nRuns

        % get functional files
        if cfg.smoothing
            funcFile = fullfile(mainPath, 'derivatives', subID, 'func', ...
                sprintf('swr%s%s_task-scenes_run-%s_bold.nii',...
                subID, 'xxxx', num2str(iRun)));
        else
            funcFile = fullfile(mainPath, 'derivatives', subID, 'func', ...
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

        % get `mcf` file
        mcf = fullfile(fmriPath, subID, 'beh', 'onsets', ...
            sprintf('mcf_%s_run-%s.mat', subID, num2str(iRun)));
        mcf_file = load(mcf);

        % get trial IDs
        trialIDs = [cell2mat(mcf_file.trialIDs),...
            ones(size(mcf_file.trialIDs))*iRun,...
            cell2mat(mcf_file.onsets)];

        % sort trial IDs and add TR
        trialIDs = sortrows(trialIDs, 3);
        trialIDs(:, 4) = round(trialIDs(:, 3) / cfg.tr) + 1;

        % remove targets
        trialIDs = trialIDs(trialIDs(:, 1) < 101, :);

        % get first and last TR of each category
        if trialIDs(1,1) < cfg.halfTrials
            firstBathroomTR = trialIDs(1, 4) - cfg.TRstartBuffer;
            lastBathroomTR = trialIDs(cfg.halfTrials, 4)...
                + cfg.TRendBuffer;
            firstKitchenTR = trialIDs(cfg.halfTrials + 1, 4)...
                - cfg.TRstartBuffer;
            lastKitchenTR = trialIDs(end, 4) + cfg.TRendBuffer;
        else
            firstKitchenTR = trialIDs(1,4) - cfg.TRstartBuffer;
            lastKitchenTR = trialIDs(cfg.halfTrials, 4)...
                + cfg.TRendBuffer;
            firstBathroomTR = trialIDs(cfg.halfTrials + 1, 4)...
                - cfg.TRstartBuffer;
            lastBathroomTR = trialIDs(end, 4) + cfg.TRendBuffer;
        end

        % get data
        v = load_untouch_nii(funcFile);
        runData = single(v.img);
        data{iRun} = runData;

        % make sure both timecourses have the same length
        timecourseLength = max([lastBathroomTR - firstBathroomTR, ...
            lastKitchenTR - firstKitchenTR]);

        % split data
        bathroomData{iRun} = runData(:, :, :, ...
            firstBathroomTR : firstBathroomTR + timecourseLength);
        kitchenData{iRun} = runData(:, :, :, ...
            firstKitchenTR : firstKitchenTR + timecourseLength);

        % check size of timecourse
        disp(['Run: ', num2str(iRun)])
        disp(['Size bathroom: ', num2str(size(bathroomData{iRun}))])
        disp(['Size kitchen: ', num2str(size(kitchenData{iRun}))])

        %% average across voxels for each mask

        bathroomMeans = {};
        kitchenMeans = {};
        for j=1:nmasks

            % get mask name
            mask_label=cfg.rois{j};
            mask_label_short = split(mask_label, '.');
            mask_label_short = mask_label_short{1};
            mask_label_short = mask_label_short(2:end);

            % check if file exists already and skip if the case
            if cfg.skipIfExists

                if cfg.smoothing
                    fileName = ['smoothed_mean_timecourse_bathroom_', mask_label_short, ...
                        '_run_', num2str(iRun), '.mat'];
                else
                    fileName = ['mean_timecourse_bathroom_', mask_label_short, ...
                        '_run_', num2str(iRun), '.mat'];
                end

                if exist(fullfile(outputdir, fileName), 'file')
                    disp(['file for ', mask_label_short, ' and ', num2str(iRun), ...
                        ' exists already'])
                    continue
                end
            end

            % check if functional or anatomical ROI
            if ismember(mask_label_short, {'PPA', 'TOS', 'RSC', 'LOC'})
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
                [1, 1, 1, size(bathroomData{iRun}, 4)]);

            % bathroom
            currentImage = bathroomData{iRun};
            filtered = currentImage(currentROI4D == 1);
            ROIImage = reshape(filtered, [], size(bathroomData{iRun}, 4));
            bathroomMeans{j} = mean(ROIImage, 1);

            % kitchen
            currentImage = kitchenData{iRun};
            filtered = currentImage(currentROI4D == 1);
            ROIImage = reshape(filtered, [], size(kitchenData{iRun}, 4));
            kitchenMeans{j} = mean(ROIImage, 1);

            % for subjects 1 in run 10 the kitchen block was not shown
            % correctly - make fMRI data nans
            if subs(iSub) == 1 && iRun == 10
                kitchenMeans{j} = nan(1, length(kitchenMeans{j}));
            end 

            % save data
            if cfg.smoothing
                fileName = ['smoothed_mean_timecourse_bathroom_', mask_label_short, ...
                    '_run_', num2str(iRun), '.mat'];
            else
                fileName = ['mean_timecourse_bathroom_', mask_label_short, ...
                    '_run_', num2str(iRun), '.mat'];
            end
            saveData = bathroomMeans{j};
            save(fullfile(outputdir, fileName), 'saveData')
            if cfg.smoothing
                fileName = ['smoothed_mean_timecourse_kitchen_', mask_label_short, ...
                    '_run_', num2str(iRun), '.mat'];
            else
                fileName = ['mean_timecourse_kitchen_', mask_label_short, ...
                    '_run_', num2str(iRun), '.mat'];
            end
            saveData = kitchenMeans{j};
            save(fullfile(outputdir, fileName), 'saveData')

        end % rois

        if cfg.saveWholeBrain && cfg.smoothing

            % save whole brain data
            fileName = ['whole_brain_timecourse_bathroom_run_', num2str(iRun)];

            % check if file exists already and skip if the case
            if cfg.skipIfExists

                if exist(fullfile(outputdir, [fileName, '.mat']), 'file')
                    disp(['file for whole brain and ', num2str(iRun), ...
                        ' exists already'])
                    continue
                end
            end

            saveData = bathroomData{iRun};
            save(fullfile(outputdir, fileName), 'saveData')
            fileName = ['whole_brain_timecourse_kitchen_run_', num2str(iRun)];
            saveData = kitchenData{iRun};
            save(fullfile(outputdir, fileName), 'saveData')
        end

    end % runs
end % subjects
end
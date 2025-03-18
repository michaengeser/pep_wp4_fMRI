function create_mcf_func(cfg)


    %% preparations
    sourcedataPath = fullfile(pwd, '..','sourcedata');
    locPath = fullfile(pwd, '..', 'localizer');
    
    % go!
    for iSub = 1:length(cfg.subNums)

        % subject ID
        subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

        %% Onsets for experimental runs
        subFolder = fullfile(sourcedataPath, subID, 'beh');
        files = dir(fullfile(subFolder, '*.tsv'));
        files = struct2table(files);
        files = sortrows(files, 'datenum','ascend');
    
        % delete sub-directory for onset files if it exists
        if exist(fullfile(subFolder, 'onsets'), 'dir')
            rmdir(fullfile(subFolder, 'onsets'),'s');
        end
    
        % make sub-directory for onset files
        mkdir(fullfile(subFolder, 'onsets'));
    
        for iRun = 1:cfg.nRuns

            fileName = [subID, '_task-main_run-', num2str(iRun), '_events.tsv'];
            runIdx = find(strcmp(fileName, files.name));
            fileData = readtable(char(fullfile(files.folder(runIdx), files.name(runIdx))),...
                'FileType', 'text', 'Delimiter', '\t');

            fileData.newTrialNum = (1:height(fileData))';

            if cfg.sortRows
                % sort table
                fileData = sortrows(fileData, 'trial','ascend');
                fileData = sortrows(fileData, 'image','ascend');
            end

            % get names, onsets and duration
            names = fileData.image;
            onsets = num2cell(fileData.Onsets);
            durations = num2cell(zeros(height(fileData), 1));
            if strcmp(fileData.category{1}, 'bathroom')
                trialIDs = num2cell(fileData.texture - 10);
            elseif strcmp(fileData.category{1}, 'kitchen')
                trialIDs = num2cell(fileData.texture - 10 + 50);
            end
            trialNums = num2cell(fileData.newTrialNum);

            % save run's onsets as `.mat` file
            fileName = sprintf('mcf_%s_run-%s.mat', subID, num2str(iRun));
            save(fullfile(subFolder, 'onsets', fileName), ...
                'names', 'onsets', 'durations', 'trialIDs', 'trialNums');
        end
    
%         %% Onsets for localizer
%         clear onsets names durations;
%         
%         % create sub-directory for onsets if it doesn't exist yet
%         if ~exist(fullfile(locPath, 'onsets'))
%             mkdir(fullfile(locPath, 'onsets'));
%         end
%     
%         % get localizer `.mat` file for current participant
%         locFile = load(fullfile(locPath, 'data', ...
%             sprintf('localizer%s.mat', num2str(cfg.subNums(iSub)))));
%     
%         % assign onsets of each image block to structure
%         i = 0;
%         catName = {'scenes' 'objects' 'scramble'};
%     
%         for cat = 1:3  % scene, object, or scramble (excl. fixation)?
%             i = i + 1;
%     
%             % get block names
%             names{i,1} = catName{cat};
%     
%             % get block onset
%             allOnsets = locFile.dat.onset( ...
%                 locFile.dat.random_trials(:,1) == cat & ...
%                 locFile.dat.random_trials(:,3) == 1 ...
%             );
%             onsets{i,1} = allOnsets;
%     
%             % get block duration
%             durations{i,1} = ones(size(onsets{i,1})) * 16;
%         end
% 
%         % save onsets
%         fileName = sprintf('mcf_%s_localizer.mat', subID);
%         save(fullfile(locPath, 'onsets', fileName), 'names', 'onsets', 'durations');
    end
end

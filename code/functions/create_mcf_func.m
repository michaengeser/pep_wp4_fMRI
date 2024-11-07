function create_mcf_func(subs)
    arguments
        subs (1,:) double
    end

    %% preparations
    subs = cellstr(string(subs));
    sourcedataPath = fullfile(pwd, '..','sourcedata');
    locPath = fullfile(pwd, '..', 'localizer');
    
    % go!
    for iSub = 1:length(subs)
        if subs(iSub) < 10
            subID = ['sub-00', num2str(subs{iSub})];
        elseif subs(iSub) < 100
            subID = ['sub-0', num2str(subs{iSub})];
        end

        %% Onsets for experimental runs
        subFolder = fullfile(sourcedataPath, subID, 'beh');
        files = dir(fullfile(subFolder, '*.tsv'));
        files = struct2table(files);
        files = sortrows(files, 'datenum','descend');
    
        % delete sub-directory for onset files if it exists
        if exist(fullfile(subFolder, 'onsets'), 'dir')
            rmdir(fullfile(subFolder, 'onsets'),'s');
        end
    
        % make sub-directory for onset files
        mkdir(fullfile(subFolder, 'onsets'));
    
        for iRun = 1:height(files)
            fileData = readtable(char(fullfile(files.folder(iRun), files.name(iRun))),...
                'FileType', 'text', 'Delimiter', '\t');

            % sort table 
            fileData = sortrows(fileData, 'texture','ascend');

            % get names, onsets and duration
            names = fileData.image;
            onsets = num2cell(fileData.trialOnset - fileData.triggerTimeStamp);
            durations = num2cell(zeros(height(fileData), 1));
            
            % save run's onsets as `.mat` file
            fileName = sprintf('mcf_%s_run-%s.mat', subID, num2str(iRun));
            save(fullfile(subFolder, 'onsets', fileName), ...
                'names', 'onsets', 'durations');
        end
    
        %% Onsets for localizer
        clear onsets names durations;
        
        % create sub-directory for onsets if it doesn't exist yet
        if ~exist(fullfile(locPath, 'onsets'))
            mkdir(fullfile(locPath, 'onsets'));
        end
    
        % get localizer `.mat` file for current participant
        locFile = load(fullfile(locPath, 'data', ...
            sprintf('localizer%s.mat', num2str(subs(iSub)))));
    
        % assign onsets of each image block to structure
        i = 0;
        catName = {'scenes' 'objects' 'scramble'};
    
        for cat = 1:3  % scene, object, or scramble (excl. fixation)?
            i = i + 1;
    
            % get block names
            names{i,1} = catName{cat};
    
            % get block onset
            allOnsets = locFile.dat.onset( ...
                locFile.dat.random_trials(:,1) == cat & ...
                locFile.dat.random_trials(:,3) == 1 ...
            );
            onsets{i,1} = allOnsets;
    
            % get block duration
            durations{i,1} = ones(size(onsets{i,1})) * 16;
        end

        % save onsets
        fileName = sprintf('mcf_%s_localizer.mat', subID);
        save(fullfile(locPath, 'onsets', fileName), 'names', 'onsets', 'durations');
    end
end

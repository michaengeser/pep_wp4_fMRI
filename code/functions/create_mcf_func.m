function create_mcf_func(subs)
    arguments
        subs (1,:) double
    end

    %% Get config.
    cfg = get_config(subs, 'create_mcf');
    
    % go!
    for iSub = 1:length(cfg.subs)
        subID = ['sub-', num2str(cfg.subs{iSub})];
    
        %% Onsets for experimental runs
        files = dir(fullfile(cfg.behPath, subID, '*.mat'));
    
        % delete sub-directory for onset files if it exists
        if exist(fullfile(cfg.behPath, subID, 'onsets'), 'dir')
            rmdir(fullfile(cfg.behPath, subID, 'onsets'),'s');
        end
    
        % make sub-directory for onset files
        mkdir(fullfile(cfg.behPath, subID, 'onsets'));
    
        for iRun = 1:length(files)
            baseFileName = files(iRun).name;
            fileData = load(fullfile(files(iRun).folder, files(iRun).name));
    
            % assign stimulus onsets to structure
            i = 0;
            for cat = 1:2  % good (= 1) or bad (= 2) image?
                for img = 1:50  % 50 images per category
                    i = i + 1;
    
                    % name images
                    names{i,1} = ['img', num2str(i)];
                    
                    % find onset of current image
                    findStim = fileData.dat.randomTrials(:,1) == cat & ...
                               fileData.dat.randomTrials(:,2) == img;
                    onsets{i,1} = fileData.dat.onset(findStim);
    
                    % code stim. duration as 0 for event-related design
                    durations{i,1} = 0;
                end
            end
            
            % save run's onsets as `.mat` file
            fileName = sprintf( ...
                'mcf_sub-%s_run-%s.mat', num2str(cfg.subs{iSub}), num2str(iRun) ...
            );
            save( ...
                fullfile(cfg.behPath, subID, 'onsets', fileName), ...
                'names', 'onsets', 'durations' ...
            );
        end
    
        %% Onsets for localizer
        clear onsets names durations;
        
        % create sub-directory for onsets if it doesn't exist yet
        if ~exist(fullfile(cfg.locPath, 'onsets'))
            mkdir(fullfile(cfg.locPath, 'onsets'));
        end
    
        % get localizer `.mat` file for current participant
        if startsWith(cfg.subs{iSub}, '0')
            fileID = extractAfter(cfg.subs{iSub}, '0');
        else
            fileID = cfg.subs{iSub};
        end
        
        locFile = load(fullfile(cfg.locPath, sprintf('localizer%s.mat', fileID)));
    
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
        fileName = sprintf('mcf_sub-%s_localizer.mat', num2str(cfg.subs{iSub}));
        save( ...
            fullfile(cfg.locPath, 'onsets', fileName), ...
            'names', 'onsets', 'durations' ...
        );
    end
end

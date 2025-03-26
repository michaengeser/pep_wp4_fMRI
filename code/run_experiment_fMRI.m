%% fMRI Experiment Script for Scene Presentation
% Written with the assistance of OpenAI's GPT-4 (2024 version).
% This script follows the BIDS standards and is designed for fMRI use with Psychtoolbox.

% Housekeeping
clear; %close all

%% Collect subject ID
cfg.subjectID = input('Enter subject number: ', 's');
cfg.runNum = input('Enter run number: ', 's');
cfg.date = datetime('today');

%% General Experiment Configuration
cfg.debugTimingFactor = 1; % Must be 1 for accurat timing (< 1 will give faster timing)
cfg.mriMode = false;  % Set to true if running in the MRI scanner, false otherwise
cfg.imageDuration = 0.25 * cfg.debugTimingFactor;  % Image presentation time in seconds
cfg.iti = 2 * cfg.debugTimingFactor;  % Inter-trial interval in seconds
cfg.startPad = 4 * cfg.debugTimingFactor;  % Time before the first trial in seconds
cfg.endPad = 10 * cfg.debugTimingFactor;  % Time after the last trial in seconds
cfg.ImageFileFormat = 'tif';

%% Paths
stimPath = fullfile(pwd,'..', 'stimuli');
outputPath = fullfile(pwd,'..', 'sourcedata', ['sub-', cfg.subjectID], 'beh');
functionPath = fullfile(pwd,'utilities');

% add functions folder to path
addpath(functionPath)

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

%% fMRI Initialization
if cfg.mriMode
    dq = InitDAQBION;  % Initialize scanner hardware (replace with actual initialization if necessary)
end

%% Image Loading
if str2double(cfg.runNum) == 0
    runCategory = 'practice';
    targetNum = str2double(cfg.runNum);
else
    if mod(str2double(cfg.runNum), 2) == 0
        runCategory = 'kitchen';
    else
        runCategory = 'bathroom';
    end
end
runImages = dir(fullfile(stimPath, runCategory, ['*.', cfg.ImageFileFormat]));
numTrials = length(runImages)*2;

%% Output File Setup
% BIDS-compliant log file
runOutputFile = fullfile(outputPath, sprintf('sub-%s_task-main_run-%s_events.tsv', cfg.subjectID, cfg.runNum));

% Check if file name exists already to avoid overwriting 
if exist(runOutputFile, 'file')
    error('FIle name exists already')
end 
fileID = fopen(runOutputFile, 'w');

% Write header for the log file
fprintf(fileID, 'subject\trun\ttrial\ttexture\tcategory\timage\ttrialOnset\titiOnset\ttrialEnd\ttriggerDate\ttriggerTimeStamp\tOnsets\tresponseKey\tresponseTime\taccuracy\n');

%% Initialize Psychtoolbox
Screen('Preference', 'SkipSyncTests', 1);  % Skip sync tests for demo purposes (remove in actual experiment)
PsychDefaultSetup(2)

% Define colors
Color.white = [255 255 255]; Color.black = [0 0 0]; Color.gray=(Color.black+Color.white)/2;
Color.red = [255 0 0]; Color.yellow= [255 255 0];

screenNumber = max(Screen('Screens'));
screencount=size(Screen('screens'),2);
if screencount>1
    windowrect=Screen(1,'rect');
    screenNumber=1; %%%%%%%%%%%%%%% 2
else
    windowrect=Screen(0,'rect');
    screenNumber=0;
end

[window, windowRect] = Screen('OpenWindow', screenNumber, [Color.gray], [], 32, 2,[], [],  kPsychNeed32BPCFloat);
[xCenter, yCenter] = RectCenter(windowRect);
ifi = Screen('GetFlipInterval', window);  % Get refresh interval
Screen('TextSize', window, 40) % define font size
DrawFormattedText(window, 'Loading...', 'center', 'center', [0 0 0]);
Screen('Flip', window);  % Show fixation cross
HideCursor; % Hide mouse cursor

% Load all images as textures
for i = 1:length(runImages)
    runTextures(i) = Screen('MakeTexture', window, imread(fullfile(runImages(i).folder, runImages(i).name)));
end

%% randomize trial order

% Ensure reproducible order across participants
rng(str2double(cfg.runNum));  % Run number as seed

% Initialize table
blkImgs = table;
blkImgs.texture = runTextures';
blkImgs.category = repmat({runCategory}, 1, length(runTextures))';
blkImgs.imgName = {runImages.name}';

% Shuffle rows
blkImgs1 = blkImgs(randperm(height(blkImgs)),:);
blkImgs2 = blkImgs(randperm(height(blkImgs)),:);
if str2double(cfg.runNum) == 0
    blkImgs = blkImgs1(1:10, :);
else
    blkImgs = [blkImgs1; blkImgs2];
end

% get targets
load(fullfile(functionPath, 'targets.mat'), 'targetStruct')
if str2double(cfg.runNum) == 0
    targetNum = numel(targetStruct);
else
    targetNum = str2double(cfg.runNum);
end

% init accuracy vector
accuracy = nan(1, numel(targetStruct(targetNum).imgName));

%% Calculate size for desired degree of visual angle

% define visual angle
x_degree = 8;
y_degree = 6;

% Get the screen resolution and viewing distance
viewing_dist = 140;
screen_hor=70; %Convert to cm (horizontal)
screen_vert=39.4; %Convert to cm (vertial)

% Calculate pixels per centimeter
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
pixPerCmX = screenXpixels / screen_hor;
pixPerCmY = screenYpixels / screen_vert;

% Calculate the size in cm for the given visual angles
sizeCmX = 2 * viewing_dist * tan(deg2rad(x_degree) / 2);
sizeCmY = 2 * viewing_dist * tan(deg2rad(y_degree) / 2);

% Convert the size from cm to pixels
sizePixX = round(sizeCmX * pixPerCmX);
sizePixY = round(sizeCmY * pixPerCmY);

% get rectangle for image of correct size
image_rect = CenterRectOnPointd([0 0 sizePixX sizePixY], xCenter, yCenter);

%% Initialize keyboard
KbName('UnifyKeyNames');
abortKey = KbName('ESCAPE');
respKeysMRI = 1:15;
% define target keys
if cfg.mriMode
    presentKey = 5;
    absentKey = 6;
else
    presentKey = 37;
    absentKey = 39;
end

%% Experiment Start
try
    % Wait for Scanner Trigger
    DrawFormattedText(window, 'Waiting for scanner...', 'center', 'center', [0 0 0]);
    Screen('Flip', window);  % Show fixation cross
    if cfg.mriMode
        [triggerDate, ~] = GetTriggerDAQBION(dq);
        triggerTimeStamp = GetSecs;
    else
        KbWait; % Wait for any key press in dummy mode
        triggerTimeStamp = GetSecs;
        triggerDate = datetime;
    end

    trialOnsets = nan(1, numTrials);  % Store trial onset times
    itiOnsets = nan(1, numTrials);  % Store ITI onset times
    trialEnd = nan(1, numTrials);  % Store trial end
    responseTimes = nan(1, numTrials);  % Store response times
    responseKeys = cell(1, numTrials);  % Store response keys
    trialAccuracy = nan(1, numTrials);  % Store response accuracy

    % Display fixation cross before the trial
    DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
    Screen('Flip', window);  % Show fixation cross

    % Wait for initial start pad (4s)
    WaitSecs(cfg.startPad);

    % Loop through all trials in the block
    for iImg = 1:height(blkImgs)

        % Initialize trial
        trialDuration = cfg.imageDuration + cfg.iti;
        responseFlag = false;
        itiFlag = false;
        responseTimes(iImg) = NaN;
        responseKeys{iImg} = 'none';
        trialAccuracy(iImg) = NaN;
        elapsedTime = 0;

        % Present image
        Screen('DrawTexture', window, blkImgs.texture(iImg), [], image_rect);
        DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
        trialOnsets(iImg) = Screen('Flip', window);  % Get trial onset time

        % trial timing
        while elapsedTime < trialDuration - ifi * 0.75

            % Check for button press (target detection task)
            [keyIsDown, responseTime, keyCode] = KbCheck;
            if cfg.mriMode
                BIONkeyCode = read(dq, 1, "OutputFormat", "Matrix");
                checkResp = BIONkeyCode(respKeysMRI);
                if sum(checkResp) >= 1
                    keyIsDown = 1;
                    responseTime = GetSecs;
                end
            end

            % Abort experiment when ESC is pressed
            if keyCode(abortKey)
                error('Experiment aborted by user')
            end

            % Show fixation cross after 250 ms
            if ~itiFlag && elapsedTime > cfg.imageDuration - ifi * 0.5
                % Inter-trial interval (ITI)
                DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Show fixation cross during ITI
                itiOnsets(iImg) = Screen('Flip', window);  % Show fixation cross
                itiFlag = true;
            end

            % Update timer
            elapsedTime = GetSecs - trialOnsets(iImg);
        end

        % Log trial end time
        trialEnd(iImg) = GetSecs;

        % check if target trial
        if ismember(iImg, targetStruct(targetNum).trialNum)

            % get target trial
            targetIdx = find(iImg == targetStruct(targetNum).trialNum);

            % check if target image is correct
            if strcmp(targetStruct(targetNum).imgName{targetIdx}, ...
                    blkImgs.imgName{iImg})
                disp('Correct target was selected')
            else
                disp(['Selected target: ', targetStruct(targetNum).imgName{targetIdx}])
                disp(['Current trial: ',  blkImgs.imgName{iImg}])
                error('Target names do not match')
            end

            % show target message
            targetMsg = ['Was the following object in the last image:', newline, ...
                newline,  targetStruct(targetNum).targetName{targetIdx}];
            DrawFormattedText(window, targetMsg, 'center', 'center', [0 0 0]);
            qTime = Screen('Flip', window);

            % Wait for response
            waitDuration = 5  * cfg.debugTimingFactor;
            responseFlag = false;
            while elapsedTime < waitDuration - ifi * 0.75
                % Check for button press (target detection task)
                [keyIsDown, responseTime, keyCode] = KbCheck;
                if cfg.mriMode
                    BIONkeyCode = read(dq, 1, "OutputFormat", "Matrix");
                    checkResp = BIONkeyCode(respKeysMRI);
                    if sum(checkResp) >= 1
                        keyIsDown = 1;
                        responseTime = GetSecs;
                    end
                end

                % Abort experiment when ESC is pressed
                if keyCode(abortKey)
                    error('Experiment aborted by user')
                end

                % Store first button press
                if keyIsDown && ~responseFlag
                    responseTimes(iImg) = responseTime;
                    if cfg.mriMode
                        responseKeys{iImg} = num2str(find(BIONkeyCode));
                    else
                        responseKeys{iImg} = num2str(find(keyCode));
                    end
                    % show fixation cross after response is provided
                    DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Show fixation cross during ITI
                    Screen('Flip', window); 
                    responseFlag = true;
                end

                % Update timer
                elapsedTime = GetSecs - qTime;
            end

            % get accuracy if response was recorded
            if strcmp(responseKeys{iImg}, 'none')
                accuracy(targetIdx) = NaN;
            else
                if targetStruct(targetNum).targetPresent(targetIdx) == 1
                    accuracy(targetIdx) = str2double(responseKeys{iImg}) == presentKey;
                else
                    accuracy(targetIdx) = str2double(responseKeys{iImg}) == absentKey;
                end
                trialAccuracy(iImg) = accuracy(targetIdx);
            end

            % Display fixation cross before the next trial
            DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
            Screen('Flip', window);  % Show fixation cross
            WaitSecs(2 * cfg.debugTimingFactor);
        end

        fprintf(fileID, '%s\t%s\t%d\t%d\t%s\t%s\t%.6f\t%.6f\t%.6f\t%s\t%.6f\t%.6f\t%s\t%.6f\t%.6f\n', ...
            cfg.subjectID, cfg.runNum, iImg,...
            blkImgs.texture(iImg), char(blkImgs.category(iImg)), char(blkImgs.imgName(iImg)), ...
            trialOnsets(iImg), itiOnsets(iImg), trialEnd(iImg), ...
            char(triggerDate), triggerTimeStamp, trialOnsets(iImg) - triggerTimeStamp,...
            responseKeys{iImg}, responseTimes(iImg), trialAccuracy(iImg)); 
    end

    % Wait for end pad
    WaitSecs(cfg.endPad);

    % show feedback message
    feedbackMsg = ['The end', newline,...
        'You answered ', num2str(round(mean(accuracy, 'omitnan')*100)), ...
        '% of the questions correctly', newline, ...
        'Thank you!'];
    DrawFormattedText(window, feedbackMsg, 'center', 'center', [0 0 0]);
    Screen('Flip', window);

    % Wait for end pad
    WaitSecs(5 * cfg.debugTimingFactor);

    % Close log file
    fclose(fileID);

    % Close Psychtoolbox Screen
    Screen('CloseAll');
    ShowCursor;

    % Show experiment duration 
    totalLength = GetSecs - triggerTimeStamp;
    disp(['The run duration was: ', char(minutes(totalLength/60))])

catch ME
    % Handle errors
    fprintf('An error occurred: %s\n', ME.message);
    Screen('CloseAll');  % Ensure the screen is closed in case of an error
    fclose(fileID);  % Close log file if open
    ShowCursor;
end

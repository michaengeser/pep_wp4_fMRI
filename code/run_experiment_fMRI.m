%% fMRI Experiment Script for Scene Presentation
% Written with the assistance of OpenAI's GPT-4 (2024 version).
% This script follows the BIDS standards and is designed for fMRI use with Psychtoolbox.

% Housekeeping
clear; close all

%% Collect subject ID
cfg.subjectID = input('Enter subject number: ', 's');
cfg.runNum = input('Enter run number: ', 's');
cfg.date = datetime('today');

%% General Experiment Configuration
cfg.mriMode = false;  % Set to true if running in the MRI scanner, false otherwise
cfg.imageDuration = 0.25;  % Image presentation time in seconds
cfg.iti = 2;  % Inter-trial interval in seconds
cfg.startPad = 4;  % Time before the first trial in seconds
cfg.endPad = 15;  % Time after the last trial in seconds
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
if mod(str2double(cfg.runNum), 2) == 0
    runCategory = 'kitchen';
else
    runCategory = 'bathroom';
end
runImages = dir(fullfile(stimPath, runCategory, ['*.', cfg.ImageFileFormat]));
numTrials = length(runImages)*2;

%% Output File Setup
% BIDS-compliant log file
runOutputFile = fullfile(outputPath, sprintf('sub-%s_task-main_run-%s_events.tsv', cfg.subjectID, cfg.runNum));
fileID = fopen(runOutputFile, 'w');

% Write header for the log file
fprintf(fileID, 'subject\trun\ttrial\ttexture\tcategory\timage\tresponseKey\tresponseTime\ttrialOnset\titiOnset\ttrialEnd\ttriggerDate\ttriggerTimeStamp\tOnsets\n');

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
blkImgs.image_name = {runImages.name}';

% Shuffle rows
blkImgs1 = blkImgs(randperm(height(blkImgs)),:);
blkImgs2 = blkImgs(randperm(height(blkImgs)),:);
blkImgs = [blkImgs1; blkImgs2];

if str2double(cfg.runNum) == 11
    blkImgs = blkImgs(1:10, :);
end 


% get targets
load(fullfile(functionPath, 'targets.mat'), 'targetStruct')

sca
error('stop')

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

    % Wait for initial start pad (4s)
    WaitSecs(cfg.startPad);

    % Display fixation cross before the trial
    DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
    Screen('Flip', window);  % Show fixation cross
    WaitSecs(2);

    % Loop through all trials in the block
    for imgNum = 1:height(currentImages)

        % Initialize trial
        trialDuration = cfg.imageDuration + cfg.iti;
        responseFlag = false;
        itiFlag = false;
        responseTimes(trialCount) = NaN;
        responseKeys{trialCount} = 'none';
        elapsedTime = 0;

        % Present image
        Screen('DrawTexture', window, currentTexture, [], image_rect);
        DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
        trialOnsets(trialCount) = Screen('Flip', window);  % Get trial onset time

        % trial timing
        while elapsedTime < trialDuration - ifi * 0.5

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
                responseTimes(trialCount) = responseTime;
                if cfg.mriMode
                    responseKeys{trialCount} = num2str(find(BIONkeyCode));
                else
                    responseKeys{trialCount} = num2str(find(keyCode));
                end
                responseFlag = true;
            end

            % Show fixation cross after 250 ms
            if ~itiFlag && elapsedTime > cfg.imageDuration - ifi * 0.5
                % Inter-trial interval (ITI)
                DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Show fixation cross during ITI
                itiOnsets(trialCount) = Screen('Flip', window);  % Show fixation cross
                itiFlag = true;
            end

            % Update timer
            elapsedTime = GetSecs - trialOnsets(trialCount);
        end

        % Log trial end time
        trialEnd(trialCount) = GetSecs;

        % Log trial information
        fprintf(fileID, '%s\t%s\t%d\t%d\t%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%.6f\n', ...
            cfg.subjectID, cfg.runNum, trialCount,...
            currentImages.texture(imgNum), char(currentImages.category(imgNum)), char(currentImages.image_name(imgNum)), ...
            responseKeys{trialCount}, responseTimes(trialCount), ...
            trialOnsets(trialCount), itiOnsets(trialCount), trialEnd(trialCount), ...
            char(triggerDate), triggerTimeStamp);

    end


    % Wait for end pad (10s)
    WaitSecs(cfg.endPad);

    % Close log file
    fclose(fileID);

    % Close Psychtoolbox Screen
    Screen('CloseAll');
    ShowCursor;

catch ME
    % Handle errors
    fprintf('An error occurred: %s\n', ME.message);
    Screen('CloseAll');  % Ensure the screen is closed in case of an error
    fclose(fileID);  % Close log file if open
    ShowCursor;
end

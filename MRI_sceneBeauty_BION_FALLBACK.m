clc; clear all; close all;

%% set paths, initialize MRI
% debug mode or actual study? Testing at home or at BION?
cfg.mriMode = 0;  % set to 1 when scanning, 0 for debugging
cfg.atBION = 0;   % set to 1 if at BION, 0 if elsewhere

if cfg.atBION == 1
    cd('D:\user\0251PF\experiment');
else
    cd('/Users/philippflieger/PhD/MRI_sceneAttractiveness/MRI_exp');
end

% add remaining paths
rng('shuffle');
addpath('./exp_stimuli');
addpath('./res');

% BION init. script
if cfg.mriMode == 1 
    dq = InitDAQBION;
end

%% define window and experimental params.
% window params.
cfg.windowColor = [0 0 0];
cfg.textColor = [255 255 255];
cfg.fixColor = [255 25 255];
cfg.fontName = 'Arial';
cfg.fontSize = 20;
[w, h] = Screen('DisplaySize', 0);  % display size in mm
cfg.screenHor = w/10;  % convert to cm
cfg.screenVert = h/10;
cfg.screenDist = 70;  % distance from screen in cm
win = Screen('Resolution', 0);
cfg.screenSize = [win.width, win.height]./2;  % get resolution

% stimulus and trial params.
cfg.picSizeHor = 600;
cfg.picSizeVert = 600;
cfg.visAngle = 7;  % `cfg.deg_object` in Daniel's code
cfg.numStims = 100;
cfg.numRuns = 8;
cfg.stimDuration = 1.45;  % `cfg.time_stim` in Daniel's code
cfg.waitPostStim = 0.1;  % time until rating; `cfg.time_wait` in Daniel's code
cfg.ITI = 0.95;
cfg.wrapWidth = w/4;  % wrap width for text presentation

% response keys for debug mode
KbName('UnifyKeyNames');
cfg.respKeys(1) = KbName('1!');
cfg.respKeys(2) = KbName('2@');
cfg.respKeys(3) = KbName('3#');
cfg.respKeys(4) = KbName('4$');

% response keys for MRI mode (CHECK THIS!!!)
cfg.respKeysMRI = [5, 6, 13, 14];

% other keys we may need
cfg.abortKey = KbName('q');
cfg.spaceKey = KbName('space');

% collect subject ID
dat.subjID = input('Enter subject number: ');
dat.run = input('Enter subject number: ');
dat.date = datetime('today');
filename = sprintf( ...
    './res/subj%s_run%s_scene-beauty_0251PF_%s.mat', ...
    num2str(dat.subjID), ...
    num2str(dat.run), ...
    dat.date ...
);

% % make sure no runs get overwritten
% while true
%     dat.run = input('Enter run number: ');
%     filename = sprintf( ...
%         './res/subj%s_run%s_res.mat', ...
%         num2str(dat.subjID), ...
%         num2str(dat.run) ...
%     );
% 
%     % check if file exists
%     if ~isfile(filename)
%         break;
%     else
%         fprintf('\n%s already exists.\n', filename);
%         overwritePrompt = input('\nOverwrite? [y/n] ', 's');
% 
%         % handle file overwriting
%         if strcmpi(overwritePrompt, 'y')
%             break;
%         else
%             correctedRunID = input('Enter CORRECT run number: ');
%             
%             if correctedRunID ~= dat.run
%                 dat.run = correctedRunID;
%                 filename = sprintf( ...
%                     './res/subj%s_run%s_res.mat', ...
%                     num2str(dat.subjID), ...
%                     num2str(dat.run) ...
%                 );
%                 break;
%             else
%                 correctedRunID = input('Enter CORRECT run number: ');
%             end
%         end
%     end
% end

%% screen setup
% if RNG doesn't work (e.g., old MATLAB version)
for i = 1:dat.subjID
    rand(1, 1);
    randperm(10);
    Shuffle(1:10);
    randi(10);
end

Screen('Preference', 'SkipSyncTests', 1);
if cfg.mriMode == 0; PsychDebugWindowConfiguration; end
[mainWindow, screenRect] = Screen('OpenWindow', 0, cfg.windowColor);
newTextSize = Screen('TextSize', mainWindow, [cfg.fontSize]);
Screen('TextFont', mainWindow, [cfg.fontName]);
cfg.monitorFrameRate = FrameRate(mainWindow);
cfg.monitorFlipInterval = Screen('GetFlipInterval', mainWindow);
sca

% adjust stim params. based on refresh rate
realStim = 1000 * cfg.stimDuration;
stimFrames = realStim / (1000/cfg.monitorFrameRate);
stimDuration = cfg.monitorFlipInterval * (stimFrames-0.5);

if cfg.mriMode ~= 0; clear hidecursor; HideCursor; end

% determine coords. of screen center
screenSize = cfg.screenSize;
center = cfg.screenSize ./ 2;
centerHor = center(1);
centerVert = center(2);

% stimulus size
picSizeCM = tan((cfg.visAngle/2) * 2 * pi/360) * cfg.screenDist;
picSizePixels = [
    picSizeCM * cfg.screenSize(1) / cfg.screenHor, ...
    picSizeCM * cfg.screenSize(1) / cfg.screenHor * ...
    (cfg.picSizeVert/cfg.picSizeHor)
];

% stimulus positioning
stimPosMid = center;
fixPos = stimPosMid;

% rectangle for image
picRect1 = [
    stimPosMid(1) - picSizePixels(1), ...
    stimPosMid(2) - picSizePixels(2)
];
picRect2 = [
    stimPosMid(1) + picSizePixels(1), ...
    stimPosMid(2) + picSizePixels(2)
];
picRect = [picRect1 picRect2];

% rectangle for stimulus
aspect = 1.5;
stimRect1 = [
    stimPosMid(1) - aspect * picSizePixels(1), ...
    stimPosMid(2) - picSizePixels(2)
];
stimRect2 = [
    stimPosMid(1) + aspect * picSizePixels(1), ...
    stimPosMid(2) + picSizePixels(2)
];
stimRect = [stimRect1 stimRect2];

% rectangle for rating wheel;
sf = 1.66;
wheelRect1 = [
    stimPosMid(1) - sf * picSizePixels(1), ...
    stimPosMid(2) - sf * picSizePixels(2)
];
wheelRect2 = [
    stimPosMid(1) + sf * picSizePixels(1), ...
    stimPosMid(2) + sf * picSizePixels(2)
];
wheelRect = [wheelRect1 wheelRect2];

%% stimulus randomization
t1 = [];
for b = 1:cfg.numRuns; t1 = [t1; randperm(cfg.numStims)']; end

for trial = 1:cfg.numStims
    stimNum = t1(trial);

    if stimNum <= cfg.numStims/2
        % good image
        dat.randomTrials(trial, 1) = 1;
        dat.randomTrials(trial, 2) = stimNum;
        dat.stim{trial} = [
            'good_image', num2str(stimNum), '.png'
        ];
    else 
        % bad image
        dat.randomTrials(trial, 1) = 2;
        dat.randomTrials(trial, 2) = stimNum - ...
            cfg.numStims/2;
        dat.stim{trial} = [
            'bad_image', num2str(stimNum - cfg.numStims/2), '.png'
        ];
    end
end

%% start presentation
% show instructions if first run
if dat.run == 0
    Screen('TextSize', mainWindow, cfg.fontSize + 10);
    headerMsg = 'Instructions\n\n\n';
    instrMsg = ...
        ['Your task is to rate images in terms of how aesthetically ' ...
        'pleasing they are to you. After an image has been ' ...
        'presented to you on the screen, select a score from ' ...
        '1 (least beautiful) to 4 (most beautiful). The scores ' ...
        'are organized around a ''rating wheel'' and their location ' ...
        'changes randomly after each image is shown (e.g., on ' ...
        'one trial, the score ''2'' is in the top left corner and on ' ...
        'the next ''2'' is in the bottom right corner). Press the ' ...
        'key on the button box corresponding to the score you''d ' ...
        'like to select.\n\n' ...
        'It is crucial that you keep your gaze fixed on the ' ...
        'center of the screen and don''t blink while an image ' ...
        'is being shown. Finally, please make sure you give your ' ...
        'rating as quickly as possible once the rating wheel is ' ...
        'on the screen. If you have any questions, please ask them ' ...
        'right now.'];

    DrawFormattedText( ...
        mainWindow, ...
        [headerMsg, instrMsg], ...  % concatenate strings
        'center', ...
        'center', ...
        cfg.textColor, ...
        cfg.wrapWidth ...
    );
    Screen('Flip', mainWindow);
    KbWait(-1);
end

% ensure correct font size is used
Screen('TextSize', mainWindow, cfg.fontSize);

% wait for scanner trigger (MRI mode) or button press (debug mode)
if cfg.mriMode == 1
    [TimeStamp, scantick] = GetTriggerDAQBION(dq);
else
    KbWait(-1); 
end

% start measuring run length
expStart = GetSecs;

% 4 seconds fix. dot before run
Screen('DrawDots', mainWindow, fixPos, 5, cfg.fixColor);
Screen('Flip', mainWindow);
WaitSecs(4);

%{
    %%%%%%%%%%%%%%%%%%%%%%
    % START PRESENTATION %
    %%%%%%%%%%%%%%%%%%%%%%
%}

ON = GetSecs;

for trial = 1:cfg.numStims
    % record time stamp at beginning of trial
    trialStartTime = GetSecs;

    % stimulus ON
    img = imread(dat.stim{trial});
    stimTex = Screen('MakeTexture', mainWindow, img);
    Screen('DrawTexture', mainWindow, stimTex, [], stimRect);
    Screen('DrawDots', mainWindow, fixPos, 5, cfg.fixColor);
    ON = Screen('Flip', mainWindow, ON + cfg.ITI);

    % stimulus OFF
    Screen('DrawDots', mainWindow, fixPos, 5, cfg.fixColor);
    ON = Screen('Flip', mainWindow, ON + cfg.stimDuration);  % `time_stim` in Daniel's code

    % draw rating wheel
    startAngle = (randi(4) - 1) * 90 + 10; % rotated 45° (prev.: `- 35;`)
    a1 = 70;
    cm = [
        20 45 80;
        90 55 168;
        165 65 165;
        255 90 255;
    ];

    for i = 1:4
        a0 = startAngle + (i-1) * a1 + (i-1) * 20;
        disp(a0);
        Screen('FillArc', mainWindow, cm(i,:), wheelRect, a0, a1);
        r = 1.2 * (wheelRect(3) - wheelRect(1)) / 2;

        % position and draw text
        xText = fixPos(1) + r * cos( ...
            deg2rad(startAngle + (i-1) * a1 + ...
            0.5 * a1 + (i-1) * 20 - 90) ...
        );
        yText = fixPos(2) + r * sin( ...
            deg2rad(startAngle + (i-1) * a1 + ...
            0.5 * a1 + (i-1) * 20 - 90) ...
        );
        DrawFormattedText( ...
            mainWindow, ...
            num2str(i), ...
            xText, ...
            yText, ...
            cfg.textColor...
        );
    end

    % wheel ON
    Screen('FillOval', mainWindow, cfg.windowColor, picRect);
    DrawFormattedText( ...
        mainWindow, ...
        'How pleasing?', ...
        'center', ...
        'center', ...
        cfg.textColor ...
    );
    ON = Screen('Flip', mainWindow, ON + cfg.waitPostStim);

    % wait for key press
    respComplete = 0;
    noAnswer = 1;
    respStartTime = GetSecs;
    tic;  % start timer

    while respComplete == 0
        if cfg.mriMode == 1
            keyCode = read(dq, 1, "OutputFormat", "Matrix");
            checkResp = keyCode(cfg.respKeysMRI);
        else
            [keyIsDown, timeSecs, keyCode] = KbCheck(-1);
            checkResp = keyCode(cfg.respKeys);
        end

        if sum(checkResp) == 1  % response in time
            dat.button(trial, 1) = find(checkResp == 1);

            % calculate RT
            respEndTime = GetSecs;
            dat.RT(trial, 1) = respEndTime - respStartTime;
            cfg.timeFeedback = 2.5 - dat.RT(trial, 1);

            % update Booleans 
            respComplete = 1;
            noAnswer = 0;
        elseif toc >= 2.5-0.1  % late response
            dat.button(trial, 1) = nan;
            dat.RT(trial, 1) = nan;
            cfg.timeFeedback = 0;
            respComplete = 1;
        elseif keyCode(cfg.abortKey) == 1  % abort experiment
            % save data and close screen
            save(filename, 'dat');
            respComplete = 1;
            sca;
        end
    end

    % response feedback
    if noAnswer == 0
        for i = 1:4
            a0 = startAngle + (i-1) * a1 + (i-1) * 20;
            Screen('FillArc', mainWindow, cm(i,:), wheelRect, a0, a1);
            r = 1.2 * (wheelRect(3) - wheelRect(1)) / 2;
    
            % position and draw text
            xText = fixPos(1) + r * cos( ...
                deg2rad(startAngle + (i-1) * a1 + ...
                0.5 * a1 + (i-1) * 20 - 90) ...
            );
            yText = fixPos(2) + r * sin( ...
                deg2rad(startAngle + (i-1) * a1 + ...
                0.5 * a1 + (i-1) * 20 - 90) ...
            );
            DrawFormattedText( ...
                mainWindow, ...
                num2str(i), ...
                xText, ...
                yText, ...
                cfg.textColor...
            );
            
            % fix wedge location after 45° rotation
            wedgeLocation = (wrapTo360(360 + a0 + a1/2) / 90) + 0.5;
    
            if wedgeLocation == dat.button(trial)
                Screen( ...
                    'FillArc', ...
                    mainWindow, ...
                    [255 255 255], ...
                    wheelRect, ...
                    a0, ...
                    a1 ...
                );
                dat.resp(trial) = i;
            elseif noAnswer == 1
                dat.resp(trial) = 999;
            end
        end
    
        Screen('FillOval', mainWindow, cfg.windowColor, picRect);
        DrawFormattedText( ...
            mainWindow, ...
            'How pleasing?', ...
            'center', ...
            'center', ...
            cfg.textColor ...
        );
        Screen('Flip', mainWindow);
    end

    % ITI
    Screen('DrawDots', mainWindow, fixPos, 5, cfg.fixColor);
    ON = Screen('Flip', mainWindow, ON + 2.5);

    % measure trial length
    trialEndTime = GetSecs;
    totalTrialTime = trialEndTime - trialStartTime;
    disp(totalTrialTime);
    dat.totalTrialTime(trial, 1) = totalTrialTime;
end

% 8 seconds fix. dot after run
Screen('DrawDots', mainWindow, fixPos, 5, cfg.fixColor);
Screen('Flip', mainWindow);
WaitSecs(8);

% get exp. length
expEnd = toc(expStart);
fprintf('Run length: %.3f seconds (%.3f minutes).', expEnd, expEnd/60);
fprintf('\nMean trial time: %.3f seconds.', mean(dat.totalTrialTime));

% save data and close screen
save(filename, 'dat');
sca;

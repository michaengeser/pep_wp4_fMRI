
%% Dynamic localizer: face, scene, object, scrambled, fixation TR 1.85
clc; clear all; close all;
%% set paths, initialize MRI
exp.atBION = 1;   % set to 1 if at BION, 0 if elsewhere
exp.mriMode = 1;

if exp.atBION == 1
    cd('D:\user\0261MK');
end

% BION init. script
if exp.mriMode == 1
    dq = InitDAQBION;
end

%% Define experiment parameters
rng('default');  % Reset MATLAB's random generator to modern settings
rng('shuffle');  % Properly seed the random generator
SubjectID = input('\nSubject ID: ','s');
time = datestr(now,'dd-mmm-yy_HH-MM');
FileName = ['DynLoc_' SubjectID '_' time];
rootDir = pwd;

%% Keys
table = KbName('KeyNames');
Key.escape =KbName('q');

%% Screen
PsychDefaultSetup(2)
% Screen('Preference', 'SkipSyncTests', 1); %skip syncronization test
% Removes the blue screen flash and minimize extraneous warnings.
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
screenNumber = max(Screen('Screens'));
% HideCursor;

% Define colors
Color.white = [255 255 255]; Color.black = [0 0 0]; Color.gray=(Color.black+Color.white)/2;
Color.red = [255 0 0]; Color.yellow= [255 255 0];

screencount=size(Screen('screens'),2);
if screencount>1
    windowrect=Screen(1,'rect');
    screenNumber=1; 
else
    windowrect=Screen(0,'rect');
    screenNumber=0;
end

[window, windowrect] = Screen('OpenWindow', screenNumber, [Color.gray], [], 32, 2,[], [],  kPsychNeed32BPCFloat); %50 50 1000 500
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
[scr.x, scr.y] = RectCenter(windowrect);
[dispWidth, dispHeight]=Screen('DisplaySize', screenNumber);

%% Stimulus

% Define block design parameters
numBlocks = 4; % Number of blocks for each type of stimuli
videosPerBlock = 6; % Number of videos to present per block
blockPerRun=8;
condMatrixShuffled = Shuffle(repmat(1:4, 1, blockPerRun));
totalVideos=60;

TR=1.85;
waitframes = 1;
blockDur = 3*videosPerBlock*numBlocks; %Each video lasts 3 sec
totalBlockdur=round(blockDur*blockPerRun);
runDur = 12*TR*numBlocks * blockPerRun + 6 * TR;  %in secs
totalDurTR= round(runDur/TR); %in TR
fprintf('\n Number of measurements = %3d TR \n ', totalDurTR);

allTrial=zeros(numBlocks*blockPerRun,videosPerBlock);
allvideoDur=zeros(numBlocks*blockPerRun,videosPerBlock);
allBlock=zeros(numBlocks*blockPerRun,1);
videoStartTimes = zeros(numBlocks*blockPerRun,videosPerBlock);

% Define the directory paths for each type of stimuli
stimuliDirs = {'face','scene','obj', 'scrambled'};
% Initialize a matrix to hold the results
selectedVideos = zeros(numBlocks, videosPerBlock);

%% Fixation task
fixationColors = [Color.red; Color.yellow]; % Define fixation colors
fixationChangeTimes = sort(randperm(totalBlockdur, 180)); % Random times for color changes
fixationIndex = 1; % To keep track of which fixation change we're on
fixationDuration = 0.5; % Duration in seconds for each color
fix.size=10;

%% Standby screen
Screen('TextFont', window, 'Arial'); Screen('TextSize', window, 50);
instText = ['Please fixate on the fixation point. \n Report the change in the color of the fixation point. \n' ...
    'Press left button when it changes to yellow. \n' ...
    'Press right button when it changes to red '];
DrawFormattedText(window, instText, 'center', 'center', Color.black,[]);
vbl=Screen('Flip', window);
ifi=Screen('GetFlipInterval', window);

% wait for scanner trigger (MRI mode) or button press (debug mode)
if exp.mriMode == 1
    [TimeStamp, scantick] = GetTriggerDAQBION(dq);
    expStart = GetSecs;
else
    KbWait(-1);
    expStart = GetSecs;
end

fixationChangeStart=expStart;
currentFixationColor=Color.black;

%% 2 TR fix before run

Screen('DrawDots', window, [scr.x scr.y], fix.size, Color.black, [], 2);
ON = Screen('Flip', window);
ON = Screen('Flip', window, ON + (2 * TR - 2*ifi));

for block = 1:numBlocks*blockPerRun %for all block
    blockStart = GetSecs;
    theCond=condMatrixShuffled(block); %choose current block
    folderPath =  fullfile(rootDir, stimuliDirs(theCond));

    for i = 1:numBlocks
        selectedVideos(i, :) = randperm(totalVideos, videosPerBlock);
    end

    for x=1:videosPerBlock
        video=selectedVideos(theCond, x);
        videoname=[stimuliDirs{theCond} '_' num2str(video) '.mp4'];
        videoPath(x) =  fullfile(folderPath, videoname);
    end

    % Open movie file
    [movie,duration, fps] = Screen('OpenMovie', window, videoPath{1});

    for trial = 1: videosPerBlock
        trialStart=GetSecs;
        % Play movie
        Screen('PlayMovie', movie, 1);
        videoStart=GetSecs;
        videoStartTimes(block, trial) = videoStart - expStart; % Normalize to fMRI scan start
        while 1
            % Fetch video frame
            tex = Screen('GetMovieImage', window, movie);

            if tex <= 0  % Check if the movie has ended
                break;
            end

            % Display the texture
            Screen('DrawTexture', window, tex);

            % Precompute fixation colors for each change time
            fixationColorsSequence = fixationColors(randi(2, length(fixationChangeTimes), 1), :);

            % Check if fixation should change
            if fixationIndex <= length(fixationChangeTimes) && GetSecs - expStart >= fixationChangeTimes(fixationIndex)
                % Set new fixation color
                currentFixationColor = fixationColorsSequence(fixationIndex, :);
                fixationChangeStart = GetSecs;
                fixationIndex = fixationIndex + 1; % Move to the next fixation time
            end

            % Determine if we should still show the new fixation color
            if GetSecs - fixationChangeStart < fixationDuration
                fixationDisplayColor = currentFixationColor; % Keep the changed color
            else
                fixationDisplayColor = Color.black; % Return to black
            end

            % Draw the fixation cross once per frame 
            Screen('DrawDots', window, [scr.x scr.y], fix.size, fixationDisplayColor, [], 2);
            ON = Screen('Flip', window);

            % Close the texture
            Screen('Close', tex);
        end

        videoEnd = GetSecs;
        videoDur = videoEnd - videoStart;
        allvideoDur(block, trial) = videoDur;

        % Stop playback
        Screen('PlayMovie', movie, 0);

        % Close movie
        Screen('CloseMovie', movie);
        nextVideoPrep=GetSecs;
        % Preload the next video if it exists
        if trial < length(videoPath)
            [movie, ~] = Screen('OpenMovie', window, videoPath{trial + 1});
        end

        trialEnd=GetSecs;
        trialdur= trialEnd-trialStart;
        allTrial(block,trial)= trialdur;
    end
    % Ensure ISI Fixation and maintain: 12 TRs
    elapsedTime = GetSecs - blockStart;  % Time spent on video
    remainingTime = (12 * TR) - elapsedTime; % Adjust ISI to fill exact 12 TRs

    Screen('DrawDots', window, [scr.x scr.y], fix.size, Color.black, [], 2);
    ON = Screen('Flip', window);
    if remainingTime > 0
        ON = Screen('Flip', window, ON + (remainingTime  - 2*ifi));
    end

    blockEnd = GetSecs;
    blockdur = blockEnd - blockStart;
    allBlock(block,1) = blockdur;
end
% 4 TR fix after run
Screen('DrawDots', window, [scr.x scr.y], fix.size, Color.black, [], 2);
ON = Screen('Flip', window);
ON = Screen('Flip', window, ON + (4 * TR - 2*ifi));

% get exp. length
expEnd = GetSecs;
expDur=expEnd-expStart;
fprintf('Run length: %.3f seconds (%.3f minutes).', expDur, expDur/60);
fprintf(' Mean Video length: %.3f seconds.',mean(mean(allvideoDur)));
fprintf('Mean Trial length: %.3f seconds.',mean(mean(allTrial)));
fprintf(' Mean Block length: %.3f seconds.',mean(allBlock));

sca;
save(char(strcat(FileName, '.mat')))
Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
KbQueueRelease();








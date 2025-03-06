%% Dynamic localizer
clc; clear all; close all;
%% set paths, initialize MRI
exp.atBION = 1;   % set to 1 if at BION, 0 if elsewhere
exp.mriMode = 1;

if exp.atBION == 1
    cd('D:\user\0250MK');
end

% BION init. script
if exp.mriMode == 1
    dq = InitDAQBION;
end

SubjectID = input('\nSubject ID: ','s');
% Seed the random number generator
rand('seed', sum(100 * clock));
time = datestr(now,'dd-mmm-yy_HH-MM');
FileName = ['DynLoc_' SubjectID '_' time];
FileName=[FileName '.txt'];
rootDir = pwd;

%% Keys
table = KbName('KeyNames');
Key.escape =27;

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
    screenNumber=1; %%%%%%%%%%%%%%% 2
else
    windowrect=Screen(0,'rect');
    screenNumber=0;
end

[window, windowrect] = Screen('OpenWindow', screenNumber, [Color.gray], [], 32, 2,[], [],  kPsychNeed32BPCFloat); %50 50 1000 500
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
[scr.x, scr.y] = RectCenter(windowrect);
[dispWidth, dispHeight]=Screen('DisplaySize', screenNumber);

scr.wid = dispWidth/10;   % cm
scr.len = dispHeight/10; % cm
scr.scr_sizeX = screenXpixels;
scr.scr_sizeY = screenYpixels;
scr.disp_sizeX = dispWidth;  %700
scr.disp_sizeY = dispHeight;  %394
scr.dist= 140; % Change distance
ifi = Screen('GetFlipInterval', window);

%% Stimulus

% Define block design parameters
numBlocks = 4; % Number of blocks for each type of stimuli
videosPerBlock = 6; % Number of videos to present per block
videoDuration = 3; % Duration of each video in seconds
blockPerRun=8;
condMatrixShuffled = Shuffle(repmat(1:4, 1, blockPerRun));
totalVideos=60;
TR=1.85;
waitframes = 1;
itiSecs=1.85*2; % blank duration in secs
itiframe=round(itiSecs/ifi); %% blank duration in frames
blockDur = videoDuration*videosPerBlock*numBlocks;
totalBlockdur=round(blockDur*blockPerRun);
meanPrep=0.5*32;
runDur=totalBlockdur+(6*TR)+meanPrep; %in secs
totalDurTR= round(runDur/TR); %in TR
fprintf('\n Number of measurements = %3d TR \n ', totalDurTR);

allTrial=zeros(numBlocks*blockPerRun,videosPerBlock);
allvideoDur=zeros(numBlocks*blockPerRun,videosPerBlock);
allBlock=zeros(numBlocks*blockPerRun,1);
allPrep=zeros(numBlocks*blockPerRun,videosPerBlock);

% Define the directory paths for each type of stimuli
stimuliDirs = {'face','scene','obj', 'scrambled'};
% stimuliDirs = {'scrambled','scrambled','scrambled', 'scrambled'};
% Initialize a matrix to hold the results
selectedVideos = zeros(numBlocks, videosPerBlock);

%% Fixation task
fixationColors = [Color.red; Color.yellow]; % Define fixation colors
fixationChangeTimes = sort(randperm(totalBlockdur, 190)); % Random times for color changes
fixationIndex = 1; % To keep track of which fixation change we're on
fixationDuration = 0.5; % Duration in seconds for each color
fix.size=pixcalc(0.25,scr);

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
else
    KbWait(-1);
end

% start measuring run length
expStart = GetSecs;
fixationChangeStart=expStart;
currentFixationColor=Color.black;

% 2 TR fix before run
for frame = 1:itiframe
    % Draw the fixation point
    Screen('DrawDots', window, [scr.x scr.y], fix.size, Color.black, [], 2);
    % Flip to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
end

KbQueueCreate; KbQueueStart;

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
    %         [movie,duration, fps] = Screen('OpenMovie', window,'D:\user\0250MK\objj\obj_1.mov');

    for trial = 1: videosPerBlock
        trialStart=GetSecs;
        % Play movie
        Screen('PlayMovie', movie, 1);
        while 1
            % Fetch video frame
            tex = Screen('GetMovieImage', window, movie);

            if tex <= 0  % Check if the movie has ended
                break;
            end

            % Display the texture
            Screen('DrawTexture', window, tex);

            % Check for fixation color change
            if fixationIndex <= length(fixationChangeTimes) && GetSecs - expStart >= fixationChangeTimes(fixationIndex)
                currentFixationColor = fixationColors(randi(2), :); % Randomly select red or yellow
                Screen('DrawDots', window, [scr.x scr.y], fix.size, currentFixationColor, [], 2);
                fixationChangeStart = GetSecs;
                fixationIndex = fixationIndex + 1; % Move to the next fixation change
            elseif GetSecs - fixationChangeStart < fixationDuration
                % Maintain the color change for the duration
                Screen('DrawDots', window, [scr.x scr.y], fix.size, currentFixationColor, [], 2);
            else
                % Draw the default fixation color outside the color change duration
                Screen('DrawDots', window, [scr.x scr.y], fix.size, Color.black, [], 2);
            end

            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

            % Close the texture
            Screen('Close', tex);
        end

        videoEnd=GetSecs;
        videoDur=videoEnd-trialStart
        allvideoDur(block,trial)= videoDur;


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
        nextVideoPrepDur=trialEnd-nextVideoPrep;
        allPrep(block,trial)=  nextVideoPrepDur;

        [pressed, fprs] = KbQueueCheck();
        if pressed && fprs(Key.escape)
            fprintf('\nKEY PRESS: Escape\n')
            Screen('CloseAll');
            KbQueueRelease;
            ShowCursor(window);
        end
    end

    blockEnd=GetSecs;
    blockdur= blockEnd-blockStart;
    allBlock(block,1)=blockdur;
end
% 4 TR fix before run
for frame = 1:(itiframe*2)
    % Draw the fixation point
    Screen('DrawDots', window, [scr.x scr.y], fix.size, Color.black, [], 2);
    % Flip to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
end

% get exp. length
expEnd = GetSecs;
expDur=expEnd-expStart;
fprintf('Run length: %.3f seconds (%.3f minutes).', expDur, expDur/60);
fprintf(' Mean Video length: %.3f seconds.',mean(mean(allvideoDur)));
fprintf('Mean Trial length: %.3f seconds.',mean(mean(allTrial)));
fprintf(' Mean Block length: %.3f seconds.',mean(allBlock));
fprintf('Mean Prep length: %.3f seconds.',mean(mean(allPrep)));


sca;
save(char(strcat(FileName, '.mat')))
Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
KbQueueRelease();








%% Define some stuff!
clc %Remove old stuff
clear all
close all

%addpath(genpath('C:\psychtoolbox\Psychtoolbox\'));

tic; %Start timer

rng('default');
rng('shuffle'); %Set random numbers to be truly random (this can throw an error -> remove it then, it does the same more clumsily below)

cfg.fmri_on=1;

cfg.windowcolor=[220,220,220];
cfg.textcolor=[50,50,50];
cfg.fixcolor=[200,50,200];
cfg.targetcolor=0.8*cfg.fixcolor;
cfg.fontname='Arial';
cfg.fontsize_instructions=26;
cfg.fontsize_fix=26;

cfg.linewidth=3;

cfg.objects={'scene_loc';'object_loc';'scrambled_loc'}; %Object (File)Names

cfg.blockorder=[randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1)]';
cfg.stim_nr=16;
cfg.trialamount=length(cfg.blockorder(:))*cfg.stim_nr;

cfg.picsize_hor=500; %Size of the images
cfg.picsize_vert=500;

cfg.deg_object=7.5; %desired degree (visual angle) of the object 

%[w,h]=Screen('DisplaySize',0); %Get display size in mm

cfg.screen_hor=70; %Convert to cm (horizontal)
cfg.screen_vert=39.4; %Convert to cm (vertial)

cfg.screen_dist=140; %Distance to screen in cm

if cfg.fmri_on==1
    w=Screen('Resolution',0); %Get screen properties
else
    w=Screen('Resolution',1); %Get screen properties
end
cfg.screensize=[w.width/2,w.height]; %Store horizontal/vertical pixels (i.e., resolution)

cfg.time_stim=0.5; %Duration of the Stimulus
cfg.time_iti=0.5; %Duration of the inter-trial interval

KbName('UnifyKeyNames'); %Fix Mac/Windows differences (hopefully)
cfg.escapeKey=KbName('escape'); %Get key Codes (indices in the KeyCode Matrix) for a few keys

%% setup
dat.subjcode=input('Enter subject number: '); %Ask for subject number

for i=1:dat.subjcode  %This is a workaround if the random number generator doesn't work (e.g., for old matlab versions) to make the random sequences different for all people
    rand(1,1);randperm(10);Shuffle(1:10);randi(10);
end

screen_ind=0;

Screen('Preference', 'SkipSyncTests', 1); %In case of problems
[mainwindow,screen_rect]=Screen('OpenWindow',screen_ind, cfg.windowcolor); %Open a fullscreen window

Screen('BlendFunction',mainwindow,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); %Initialize alpha blending (needed for making stuff transparent, see part 2 of the experiment)

Screen('TextSize', mainwindow, [cfg.fontsize_instructions]);  %Specify Font Size for instructions
Screen('TextFont', mainwindow, [cfg.fontname]);  %Specify Font

cfg.monitorFrameRate=FrameRate(mainwindow); %Get the frame rate of the screen
cfg.monitorFlipInterval=Screen('GetFlipInterval',mainwindow); %Get the interval between flips (is a function of the framerate, e.g., 16ms for 60hz)

%Here, we convert the timing-crucial frames into a number of flips. That
%is, we round the timing to the closest flip (e.g., a time of 36ms cannot
%be reached with a 60hz screen, but only 16ms, 33ms, 50ms, etc. Because of
%that we need to make sure we get the best we can (which in this example is
%33. This is what these lines do.
real_stim=1000*cfg.time_stim; %For the stimulus
cfg.stim_frames=real_stim/(1000/cfg.monitorFrameRate);
time_stim=cfg.monitorFlipInterval*(cfg.stim_frames-0.5);

real_stim=1000*cfg.time_iti; %For the iti
cfg.iti_frames=real_stim/(1000/cfg.monitorFrameRate);
time_iti=cfg.monitorFlipInterval*(cfg.iti_frames-0.5);

HideCursor;

%% initialize MRI stuff
if cfg.fmri_on==1 
    dq=InitDAQBION;
end

%% rectangles
%Now we convert the desired degree visual angle into pixels. First we
%compute how many cm the image should subtend at a certain visual angle.
%Then we convert the cm size into pixel size
picsizecm=tan((cfg.deg_object)*2*pi/360)*cfg.screen_dist;  %Centimeters size of the picture
pic_size=[picsizecm*cfg.screensize(1)/cfg.screen_hor,picsizecm*cfg.screensize(1)/cfg.screen_hor];  %Pixels of the picture

%For the (center) location rectangle
center = cfg.screensize./2;  %Determine center coordinates (center)
picrectangle1=[center(1)-pic_size(1)/2,center(2)-(pic_size(2)/cfg.picsize_hor*cfg.picsize_vert)/2]; %Top left corner
picrectangle2=[center(1)+pic_size(1)/2,center(2)+(pic_size(2)/cfg.picsize_hor*cfg.picsize_vert)/2]; %Bottom right corner
cfg.rectangle=[picrectangle1,picrectangle2]; %The two corners in the format [x1,y1,x2,y2] make the rectangle

%fix rectangle
picrectangle1=[center(1)-pic_size(1)/2/12,center(2)-(pic_size(2)/cfg.picsize_hor*cfg.picsize_vert)/2/12];
picrectangle2=[center(1)+pic_size(1)/2/12,center(2)+(pic_size(2)/cfg.picsize_hor*cfg.picsize_vert)/2/12];
cfg.rectangle_fix=[picrectangle1,picrectangle2];

%% randomization
cfg.blockorder=cfg.blockorder(:);
while sum(diff(cfg.blockorder)==0)>0
    cfg.blockorder=[randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1); randperm(length(cfg.objects)+1)]';
    cfg.blockorder=cfg.blockorder(:);
end
dat.random_trials=[];
for b=1:length(cfg.blockorder)
    current_block(:,1)=repmat(cfg.blockorder(b),cfg.stim_nr,1);
    current_block(:,2)=randperm(cfg.stim_nr);
    current_block(:,3)=[1;zeros(cfg.stim_nr-1,1)];
    current_block(:,4)=[Shuffle([1;zeros(cfg.stim_nr/2-1,1)]);Shuffle([1;zeros(cfg.stim_nr/2-1,1)])];
    dat.random_trials=[dat.random_trials;current_block];
end

%% set positions for fixation cross
xy_fix=[];
xy_fix(1,:)=[cfg.rectangle_fix(1),cfg.rectangle_fix(3) ...
             cfg.screensize(1)/2,cfg.screensize(1)/2];
xy_fix(2,:)=[cfg.screensize(2)/2,cfg.screensize(2)/2, ...
             cfg.rectangle_fix(2),cfg.rectangle_fix(4)];
         
%% Wait for scanner
DrawFormattedText(mainwindow,'Waiting for scanner trigger...','center','center',cfg.textcolor);
Screen('Flip',mainwindow);

if cfg.fmri_on==1 %put scanner trigger lines here
    [TimeStamp,scantick]=GetTriggerDAQBION(dq);
else
    KbWait(-1); 
end
fmri_start=GetSecs;

%% trial loop
for trial=1:cfg.trialamount %The actual trial loop where stuff happens

    if trial==1
        %draw grid & fixation
        Screen('FillRect',mainwindow,[255,255,255],cfg.rectangle)
        Screen('DrawLines',mainwindow,xy_fix,cfg.linewidth,cfg.fixcolor);  
        Screen('FrameRect',mainwindow,cfg.textcolor,cfg.rectangle,cfg.linewidth)
        ON=Screen('Flip',mainwindow);
    end

    %draw stimulus
    if dat.random_trials(trial,1)<max(cfg.blockorder) 
        dat.stim{trial}=['./stimuli/',cfg.objects{dat.random_trials(trial,1)},num2str(dat.random_trials(trial,2))];
        stim_tex=Screen('MakeTexture',mainwindow,imread([dat.stim{trial},'.png'])); %Create Texture
        Screen('DrawTexture',mainwindow,stim_tex,[],cfg.rectangle); %Draw Stimulus
    else
        Screen('FillRect',mainwindow,[255,255,255],cfg.rectangle)
    end 

    %draw grid & fixation
    Screen('DrawLines',mainwindow,xy_fix,cfg.linewidth,cfg.fixcolor); 
    Screen('FrameRect',mainwindow,cfg.textcolor,cfg.rectangle,cfg.linewidth);
        
    %flip
    if trial>1
        ON=Screen('Flip',mainwindow,ON+time_iti);
    else
        ON=Screen('Flip',mainwindow,fmri_start+4);
    end 
    dat.onset(trial)=ON-fmri_start; 

    %check for keyboard
    KeyIsDown=0;
    [keyIsDown,timeSecs,keyCode]=KbCheck(-1); 

    %code response
    if keyIsDown==1
        dat.resp(trial)=1;
    end

    %If you pressed the escape key, it will exit the loop here
    if keyCode(cfg.escapeKey)==1 
       sca 
    end
    
    %blank the stimulus
    Screen('FillRect',mainwindow,[255,255,255],cfg.rectangle)
    Screen('DrawLines',mainwindow,xy_fix,cfg.linewidth,cfg.fixcolor);  
    Screen('FrameRect',mainwindow,cfg.textcolor,cfg.rectangle,cfg.linewidth)
    ON=Screen('Flip',mainwindow,ON+time_stim);
  
    if dat.random_trials(trial,1)<max(cfg.blockorder)
        Screen('Close',stim_tex); %Close textures
    end
    
end

WaitSecs(10-time_iti);
sca %Close the window
dat.total_time=GetSecs-fmri_start; %Get total time

%% Save
save(['./data/localizer',num2str(dat.subjcode)],'cfg','dat'); %Save the data

%
function Drift_Loom
clearvars -global -except Ard; clearvars; sca;
% clear all; clear global; clear mex; sca; %#ok<CLMEX,CLALL>
[scriptPath, scriptName] = fileparts(mfilename('fullpath')); %#ok<*ASGLU>
TimeAnimation=0; stim_id_list=[];
Screen('Preference', 'SkipSyncTests', 1); %to skip synchronization failure set to 1

calc_frames = 1 ; frame_time = 0.114;%0.113%0.0088%0.044%0.176%0.264;
%define target eye (1 is left eye or single monitor/dichoptic=0. 2 is right eye)
tar_eye = 2;
Param.tar_eye = tar_eye;
Param. Dichoptic = 1 ; Param.stereoMode = 0; %1
Param. screenid  = 0 ; % it matters only if Param.Dichoptic=0;
monitor_dist_cm  = 13;%13;
rot_bias_left   = -0 ; % if eye_rot_bias>0, set [+x,-x] L,R; else if eye_rot_bias<0, set [-x,+x] L,R
% rot_bias_right  = +0 ;
% Set random or sequential sequence of stimuli (directions or positions):
Param.seqmode='random'; %'sequential'
% Param.seqmode='sequential';
Param.fullyInterleaved = 1 ; %
Param.prestimInterval  = 30 ; %sec
Param.usePhotodiode = 1;
Param.useDAQdevice  = 0;
Param.useArduino    = 1;
test = 0 ;          % set to 1 if you want to test the stimulation without external triggering
test_screen = 0 ;   % set to 1 if you want to test the stimulation on a smaller screen; 0=full field.
loadMyGammaTable = 0 ;

Param.StimProtocols.   DriftingGrating           = 2 ; % set to 2 for "navigation stim"
Param.StimProtocols.   Looming                   = 0 ;
Param.StimProtocols.   RFMapping                 = 0 ;



%% DG
reps       = [10];
stim_time       = [5];
tf         = 2;%[1 2 4];
sf         = 0.04;%[0.02 0.04 0.08];
directions = 12;
angles     = [];%[0 90 180 270];
if Param.StimProtocols.DriftingGrating
    global DG %#ok<*TLEV>
    DG.n_reps = reps(1) ; % number of repetitions (1 repetition = each direction once).
    DG.drawMask = 0; % draw grating through circular aperture
    DG.maskCenter_offset = [+0/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
    DG.maskSize_deg = 70; %deg
    DG.spacfreq = sf; %cpd
    DG.cyclespersecond = tf;
    DG.directions = directions(1); % Set the number of directions you want to be presented:
    DG.angles = angles;%[45 90 225 270];%[]; % specify cartesian angles
    DG.offset_rot_deg = [  rot_bias_left  ] ; %#ok<*NBRAK>
    DG.stimulus_time = stim_time(1) ; % seconds per stimulus direction
    DG.poststim_time = 3 ; % seconds per blank stimulus (gray screen that precedes the stimulus)
    DG.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG.phase = []; % initial phase
    DG.amplitude = 0.30; % contrast (if 0.5->100%)
    DG.BackgroundLuminance = 127;
    DG.stim_id = 1;  stim_id_list=[stim_id_list;DG.stim_id];   % for the 'putsample' value
    if Param.StimProtocols.DriftingGrating==2
        DG.n_reps = 300 ;
        DG.drawMask = 0;
        DG.spacfreq = 0.04; %cpd
        DG.cyclespersecond = 2 ;
        DG.directions = 12; DG.angles = [];
        DG.stimulus_time = 2 ;
        DG.poststim_time = 3 ;
        Param.prestimInterval  = 1 ; %sec
    end
    DG.time = DG.n_reps*DG.directions*length(DG.spacfreq)...
        *length(DG.cyclespersecond)*(DG.poststim_time+DG.stimulus_time);
    TimeAnimation = TimeAnimation + DG.time;
end
if Param.StimProtocols.Looming
    global LO %#ok<*TLEV>
    LO.n_reps = 10; % number of repetitions (1 repetition = each direction once).
   
%     LO.stimulus_time = stim(1) ; % seconds per stimulus direction
    LO.poststim_time = 3 ; % seconds per blank stimulus (gray screen that precedes the stimulus)
%     LO.expansion_speed = 1.03;%factor of enlargement of the circle
%     LO.positions = 8;%number of positions to sample
    LO.expansion_speeds = [5 10 20 40 80 160];%in degrees per second
    LO.colors = [0 0 0;255 255 255];
    LO.BackgroundLuminance = 127;
    LO.stim_id = 1.3;  stim_id_list=[stim_id_list;LO.stim_id];   % for the 'putsample' value
    if Param.StimProtocols.Looming==2
        LO.n_reps = 300 ;
%         LO.stimulus_time = 2 ;
        LO.poststim_time = 3 ;
        Param.prestimInterval  = 1 ; %sec
    end
%     LO.time = LO.n_reps*(LO.poststim_time+LO.stimulus_time);
%     TimeAnimation = TimeAnimation + LO.time;
end
if ismember(1,Param.StimProtocols.RFMapping)
    global RFM
    RFM.n_reps = 10 ;
    RFM.patch_time = 0.5;
    RFM.interpatch_time = 0.5;%0.030 ; % seconds bewteen patches.
    RFM.n_patches = [10, 10];
    RFM.fov = [80,80];  % to obtain patches of 10deg %[96,64];   %[96,80]; [84,72]  % size of the field of view in degree
    RFM.view_offset = [0 0];%[-17,-19];        % offset of the field of view (x<0 => sx , y<0 => up)
    RFM.rel_patch_size = 1.2;%1.2;           % patch size: 1: touching  - 0.5: size an distance is equal  %Ale: relative size, if 1 tha patches touch precisely, >1 they overlap
    RFM.Black = 1; RFM.BlackLuminance =  0;
    RFM.White = 1; RFM.WhiteLuminance = 255;
                   RFM.BackgroundLuminance = 127;
    RFM.gnomonic=0;
    RFM.stim_id = 3.6;  stim_id_list=[stim_id_list;RFM.stim_id];
    RFM.time = RFM.n_reps*prod(RFM.n_patches)*2*(RFM.patch_time+RFM.interpatch_time);
    TimeAnimation = TimeAnimation + RFM.time;
%     if Param.StimProtocols.DriftingGrating || Param.StimProtocols.RetinotopicPosition || Param.StimProtocols.PhiGratings
%         error('RFMapping must be the only StimProtocol ');
%     end
    RFM.save_screen  = 0 ;
    show_screen_only = 0 ; % if 1 save and show the screen masks and exit without presenting any stimuli.
end
%% Set screen parameters, saving folder, DAQ device
% Choose screen:
% screenid = max(Screen('Screens'));
numScreens = Screen('Screens');
if Param.Dichoptic
    if Param.stereoMode == 4 || Param.stereoMode == 1
        if numel(numScreens) == 3
            Param.screenid(1) = 1 ;
        else
            Param.screenid(1) = 1 ;
        end
            mouse_dist_cm([1 2]) = [monitor_dist_cm monitor_dist_cm] ;
            screenSize_cm_horz = 53;%44;
    elseif Param.stereoMode == 0
%         Param.screenid(1) = 1 ;
        Param.screenid(1) = 2 ;%DRAGO
        Param.screenid(2) = 3 ;%DRAGO
        mouse_dist_cm([1 2]) = [monitor_dist_cm monitor_dist_cm] ;
        screenSize_cm_horz = 44;
      
%         Param.screenid(1) = 1 ;
%         Param.screenid(2) = 2 ;
%         mouse_dist_cm(1) = 21 ;
%         mouse_dist_cm(2) = 18 ;
%         screenSize_cm_horz = 44;
    end

else %Param.Dichoptic==0
    if Param.screenid == 1
        mouse_dist_cm = 13;%20%13.3%#ok<*NOPRT> %12.7%19 %18.8 ;
        screenSize_cm_horz = 44;%47.5%44;
    elseif Param.screenid == 2
        mouse_dist_cm = 18;%20 %18.8 ;
        screenSize_cm_horz = 44;
    elseif Param.screenid == 0
        %         mouse_dist_cm = 20%18.8 %18.8 ;
        %         screenSize_cm_horz = 88;
        screenSize_cm_horz = 44;
    end
end
%% Open window and get basic screen parameters

% Open a fullscreen onscreen window on that display:
% if test_screen==0 && ( Param.StimProtocols.RandomDotStereoEdge || Param.StimProtocols.RandomDotStereoBar )
%     oldRes = Screen('Resolution', Param.screenid(1));
%     if oldRes.width ~= 1600 || oldRes.height ~= 900
%         disp(' '); disp(' ');
%         disp('  Warning:  Changing screen resolution to 1600x900'); pause(2);
%         oldRes = SetResolution(Param.screenid(1), 1600, 900);
%     end
% end
if     test_screen == 0
    if loadMyGammaTable
        % get calibrated gammaTable variable:
        if     strfind(eval('computer'),'WIN')
            load('C:\Users\lachioma\Dropbox\Code\alessandro\screen_calibration\2016-05-04\MyGammaTable.mat');
        elseif strfind(eval('computer'),'LNX')
            load('/home/alessandro/Dropbox/Code/alessandro/screen_calibration/2016-05-04/MyGammaTable.mat');
        end
    end
    if Param.Dichoptic
        if Param.stereoMode == 10
            [screen(1), Param.screenRect ] = Screen('OpenWindow', Param.screenid(1), 200, [],[],[], Param.stereoMode);
            slaveScreen = 2;
            [screen(2), Param.screenRectR] = Screen('OpenWindow', slaveScreen, 200, [], [], [], Param.stereoMode);
            if loadMyGammaTable
                %                 Screen('LoadNormalizedGammaTable', screen(1), gammaTable1*[1 1 1]);
                %                 Screen('LoadNormalizedGammaTable', screen(2), gammaTable2*[1 1 1]);
            end
        elseif Param.stereoMode == 4 || Param.stereoMode == 1
            [screen(1), Param.screenRect ] = PsychImaging('OpenWindow', Param.screenid(1), 200, [],[],[], Param.stereoMode);
            if loadMyGammaTable
                Screen('LoadNormalizedGammaTable', screen(1), gammaTable*[1 1 1]);
            end
        elseif Param.stereoMode == 0
%             [screen(1), Param.screenRect ] = PsychImaging('OpenWindow', Param.screenid(1), 200, [],[],[], Param.stereoMode);
            %             [screen(1), Param.screenRect ] = Screen('OpenWindow', Param.screenid(1), 200 );
            %             [screen(2), Param.screenRectR] = Screen('OpenWindow', Param.screenid(2), 200 );
            [screen(1), Param.screenRect ] = Screen('OpenWindow', Param.screenid(1), 200 ); %DRAGO
            [screen(2), Param.screenRectR] = Screen('OpenWindow', Param.screenid(2), 200 ); %DRAGO
       
            if loadMyGammaTable
                Screen('LoadNormalizedGammaTable', screen(1), gammaTable1*[1 1 1]);
                Screen('LoadNormalizedGammaTable', screen(2), gammaTable2*[1 1 1]);
            end
        end
    else % Dichoptic==0
        [screen(1), Param.screenRect] = Screen('OpenWindow', Param.screenid(1), 200 );
        if loadMyGammaTable
            if     Param.screenid(1) == 1
                Screen('LoadNormalizedGammaTable', screen(1), gammaTable1*[1 1 1]);
            elseif Param.screenid(1) == 2
                Screen('LoadNormalizedGammaTable', screen(1), gammaTable2*[1 1 1]);
            end
        end
    end
    HideCursor;
elseif test_screen == 1
    mon_pos = get(0,'MonitorPositions');
    pos_screen1 = [         50       50  round(mon_pos(1,3)*1/2)                  round(mon_pos(1,4)*2/3) ];
    %     pos_screen2 = [ mon_pos(1,1)+50  50  mon_pos(1,1)+50+round(mon_pos(1,3)*2/3)  round(mon_pos(1,4)*2/3) ];
    %     if Param.Dichoptic
    %         if Param.stereoMode == 4 || Param.stereoMode == 1
    %             [screen(1), Param.screenRect ] = PsychImaging('OpenWindow', Param.screenid(1), 200, pos_screen1,[],[], Param.stereoMode);
    % %             [screen(1), Param.screenRect ] = Screen('OpenWindow', Param.screenid(1), 127, pos_screen1,[],[], Param.stereoMode);
    %         elseif Param.stereoMode == 0
    %             [screen(1), Param.screenRect ] = Screen('OpenWindow', Param.screenid(1), 127, pos_screen1 );
    %             [screen(2), Param.screenRectR] = Screen('OpenWindow', Param.screenid(2), 127, pos_screen2 );
    %         end
    %     else %Param.Dichoptic==0
    [screen(1), Param.screenRect] = Screen('OpenWindow', Param.screenid(1), 127, pos_screen1 );
    %     end
end
if loadMyGammaTable == 0
%     LoadIdentityClut(screen);
    for screen_var = 1:length(screen)

        LoadIdentityClut(screen(screen_var));
    end
end
% when stereoMode==4 or 1, Param.screenRect(3) is half the whole screen width
Param.screenRes = [Param.screenRect(3), Param.screenRect(4)];


try
    % If Param.useDAQdevice==1 set DAQ device for stim analog output and digital input from laser shutter:
    SetDAQdevice_64bit
catch
    Screen('CloseAll');
    Priority(0);
    ShowCursor;
    % Bye bye!
    clearvars -global -except Ard
    clearvars;
    return;
end
if Param.usePhotodiode
    Param.Photodiode = SetPhotodiode(Param.screenRect);
    % winPD = screen(1);
    % Screen('FillRect', winPD , Param.Photodiode.ColorOff, Param.Photodiode.Coord );
end

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;
% Make sure the GLSL shading language is supported:
AssertGLSL;
Priority(MaxPriority(screen));
% Retrieve video redraw interval for later control of our animation timing:
Param.ifi = Screen('GetFlipInterval', screen(1));
% Calculates the visual angle subtended by a single pixel:
% see http://psychtoolbox.org/PTB-2/download/FineTutorial/revision_2/MatClassAll.rtf
Param.pixel_cm  = screenSize_cm_horz/Param.screenRect(3); % calculates the size of a pixel in cm
degperpix = (2*atan(Param.pixel_cm./(2*mouse_dist_cm(1)))).*(180/pi);  %%true only for small angles, more generally would be: atan(pix./animal_dist_cm)  %*(180/pi) is necessary because atan in matlab gives result in radians
Param.pixperdeg = 1./degperpix;
% if Param.Dichoptic
%     degperpix = (2*atan(Param.pixel_cm./(2*mouse_dist_cm(2)))).*(180/pi);  %%true only for small angles, more generally would be: atan(pix./animal_dist_cm)  %*(180/pi) is necessary because atan in matlab gives result in radians
%     Param.pixperdeg(2) = 1./degperpix;
% end
Param.mouse_dist_cm = mouse_dist_cm;
Param.screenSize_cm_horz = screenSize_cm_horz;
%% Count frames and set total times 

%if the looming is activated
if Param.StimProtocols.Looming
    %calculate the stimulus time based on the size of the screen and the
    %specified angular expansion speeds
    %get 9100% of the screen y size in degrees
    screen_angularSize = Param.screenRect(4)*degperpix;
    
    %calculate times of the stimulation based on the expansion speed
    LO.time_perstimulus = screen_angularSize./LO.expansion_speeds;
    LO.time = sum(LO.time_perstimulus + LO.poststim_time).*LO.n_reps.*size(LO.colors,1);
    TimeAnimation = TimeAnimation + LO.time;
end

TimeAnimation = TimeAnimation + Param.prestimInterval;
framesTOT = ceil( TimeAnimation / frame_time );
disp(['Total time for whole stimulation is ' num2str(TimeAnimation/60) ' min.'])
disp(['Total number of frames with frame time ' num2str(frame_time) ' is ' num2str(framesTOT) ])
if calc_frames
    Screen('CloseAll');
    Priority(0);
    ShowCursor;
    return
end
% Check if the stim_id of the stimuli are far enough:
if Param.useDAQdevice==1
    thr_stimid=0.29;
else
    thr_stimid=0;
end
if any(abs(diff(sort(stim_id_list))) <= thr_stimid )
    warndlg({['The stim_id used are too close, you won''t be able to distinguish the different stimuli! Aborting...'];...
        [ 'stim_id_list = ' mat2str(stim_id_list')]} )
    return
end

if     strfind(eval('computer'),'WIN')
    save_dir=['C:\Users\drguggiana\Dropbox\Stimuli\Stim_data\' datestr(now, 'yyyy_mm_dd')];
elseif strfind(eval('computer'),'LNX')
    save_dir=['/home/alessandro/Dropbox/Code/alessandro/stim_data/' datestr(now, 'yyyy_mm_dd')];
end
if ~isdir(save_dir)
    mkdir(save_dir)
end
datetime_suffix = datestr(now, 'yyyy_mm_dd_HH_MM_SS') ;
global stim_fname
stim_fname = [save_dir filesep 'InfoStim-' datetime_suffix '-DriftGrating_RetPos_PhiMotion'];
Param.save_dir = save_dir; Param.datetime_suffix = datetime_suffix;
Param.stim_fname = stim_fname;

internalRotation = 1; % 1=rotation of the gratings inside the "box", 0=rotation of the box.
if   internalRotation
    Param.rotateMode = kPsychUseTextureMatrixForRotation;
else
    Param.rotateMode = [];
end
%% Set parameters of stimulations
Param.stimSeq = {};
if Param.StimProtocols.DriftingGrating
    %     if Param.StimProtocols.RetinotopicPosition
    %         DG.BackgroundLuminance = RP.BackgroundLuminance ;
    %     end
    [ DG, Param ] = SetDriftingGratings( screen(tar_eye), DG, Param )
end
if Param.StimProtocols.Looming
    %     if Param.StimProtocols.RetinotopicPosition
    %         DG.BackgroundLuminance = RP.BackgroundLuminance ;
    %     end
    [ LO, Param ] = SetLooming( screen(tar_eye), LO, Param )
end
if ismember(1, Param.StimProtocols.RFMapping)
%     win=screen(1);
    [ RFM, Param ] = SetRFMapping( screen(tar_eye), RFM, Param )

%     if show_screen_only==1
%         for i=1:length(screen)
%             mon_pos = get(0,'MonitorPositions');
%             figure; set(gcf,'Position',[mon_pos(2,1),mon_pos(2,2),500,500]);
%             imshow(RFM.ScreenImg_merge(:,:,:,i),'InitialMagnification','fit','Border','tight'); %#ok<*NODEF>
%         end
%         Screen('CloseAll'); Priority(0); ShowCursor;
%         return
% %     else commandwindow;
%     end
end
% make stimSeq a row numeric array:
Param.stimSeq = cat(2, Param.stimSeq{:} );
if Param.fullyInterleaved
    %if there is RF mapping
    if ismember(1, Param.StimProtocols.RFMapping)
        %interleave only the non RF mapping trials. put RF mapping at the end
        rf_trials = Param.stimSeq(Param.stimSeq==RFM.stim_id);
        other_trials = Param.stimSeq(Param.stimSeq~=RFM.stim_id);
        Param.stimSeq = cat(2,other_trials(randperm(length(other_trials))),rf_trials);
    else
        Param.stimSeq = Param.stimSeq(randperm(length(Param.stimSeq))); % elements 1 to 4 shuffled
    end
end

Param.test=test;
Param.test_screen=test_screen;
Param.prestimInterval_fr=round(Param.prestimInterval/Param.ifi);
save(stim_fname, 'Param');

%% Present stimuli
prestimColor = 127;
if Param.StimProtocols.DriftingGrating
    prestimColor = DG.BackgroundLuminance ;
end
if Param.StimProtocols.Looming
    prestimColor = LO.BackgroundLuminance ;
end
if Param.StimProtocols.RFMapping
    prestimColor = RFM.BackgroundLuminance ;
end

if     Param.stereoMode==0
    for i=1:length(screen)
        Screen('FillRect', screen(i), prestimColor);
        Screen('Flip', screen(i));
    end
    % elseif Param.stereoMode==4 || Param.stereoMode == 1
    %     Screen('SelectStereoDrawBuffer', screen(1), 0);
    %     Screen('FillRect', screen(1), prestimColor);
    %     Screen('SelectStereoDrawBuffer', screen(1), 1);
    %     Screen('FillRect', screen(1), prestimColor);
    %     Screen('Flip', screen(1));
end
% Wait trigger to start stimulation protocol:
if test == 1 || Param.useDAQdevice==0 % Wait for release of all keys on keyboard, then sync us to retrace:
    %     startKey = KbName('space'); keyCode = [];
    % 	while ~keyCode(startKey)  %Ale: keep pausing as long as you don't press the startKey
    %        [~,~,keyCode] = KbCheck;
    %        KbReleaseWait;
    disp(['Total time for whole stimulation is ' num2str(TimeAnimation/60) ' min.'])
    disp(['Total number of frames with frame time ' num2str(frame_time) ' is ' num2str(framesTOT) ])
    disp(' Waiting for starting input..... Press Ctrl to start!');
    while ~KbCheck
    end
    pause(0.5)
else  % Wait for trigger
    disp(' '); disp(' ');
    disp(' Waiting for trigger input.....');
    while inputSingleScan(nidaq)<0.5;
        Exit_If_Shift(test,nidaq,'only_shift')
    end
end
disp(' ');
disp([' ---> Animation started!  ( preStim interval ' num2str(Param.prestimInterval) 's)']);
pause(1); disp(' ');
% Initial prestim interval with grey screen
% just a flip to get the vbl:
vbl = Screen('Flip', screen(1));

tic_start=tic;
for i = 1:Param.prestimInterval_fr
    if     Param.stereoMode==0
        for i=1:length(screen)
            vbl = Screen('Flip', screen(i), vbl + 0.5 * Param.ifi);
        end
    else
        vbl = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi);
    end
    Exit_If_Shift(test,nidaq)
end
% Animation loop: press 'shift' to exit
% (stim files are saved when you interrupt)
for s = 1:length(Param.stimSeq) % for all stimuli and reps
    %     tic_stim=tic;
    if 	   exist('DG','var') && Param.stimSeq(s)==DG.stim_id %
        [DG,  vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, DG , tar_eye);
        DG.cnt = DG.cnt+1;    
    elseif exist('LO','var') && Param.stimSeq(s)==LO.stim_id %
        [LO,  vbl] = PresentLooming( nidaq, screen, vbl, Param, LO , tar_eye);
        LO.cnt = LO.cnt+1;
    elseif exist('RFM','var') && Param.stimSeq(s)==RFM.stim_id %
        [RFM, vbl] = PresentRFMapping( nidaq, screen, vbl, Param, RFM , tar_eye);
        RFM.cnt = RFM.cnt+1;
    end
    if exist('scriptName','var') && s==1
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
end

StimSignal( Param, test, nidaq, 0 )

if     Param.stereoMode==0
    for i=1:length(screen)
        Screen('FillRect', screen(i), 80);
        Screen('Flip', screen(i));
    end
    % elseif Param.stereoMode==4 || Param.stereoMode == 1
    %     Screen('SelectStereoDrawBuffer', screen(1), 0);
    %     Screen('FillRect', screen(1), 80);
    %     Screen('SelectStereoDrawBuffer', screen(1), 1);
    %     Screen('FillRect', screen(1), 80);
    %     Screen('Flip', screen(1));
end

%% Closes windows and save parameters:

% SaveParameters_script;
SaveSimple_3;
disp(['* StimInfo file saved in ' save_dir ]);


Param.tempo_tot=toc(tic_start);
disp(['* Whole animation lasted ' num2str(Param.tempo_tot/60) ' min.' ]);
if test == 1 || Param.useDAQdevice==0 % Wait key press to close screen:
    while KbCheck == 0
    end
elseif  test == 0 && Param.useDAQdevice==1 % Wait shutter closure to close screen:
    while inputSingleScan(nidaq)>0.5
    end
end
LoadIdentityClut(screen);
Screen('CloseAll');
Priority(0);
ShowCursor;


% % catch ME
% % 	% ---------- Error Handling ----------
% % 	% If there is an error in our code, we will end up here.
% %
% % 	% The try-catch block ensures that Screen will restore the display and return us
% % 	% to the MATLAB prompt even if there is an error in our code.  Without this try-catch
% % 	% block, Screen could still have control of the display when MATLAB throws an error, in
% % 	% which case the user will not see the MATLAB prompt.
% % 	% We throw the error again so the user sees the error description.
% %     if strcmp(ME.identifier,'MATLAB:dispatcher:InexactCaseMatch') &&...
% %         ~isempty(strfind(ME.stack(1).name,'Exit_If_Shift'))
% %         % nothing to do
% %         disp('   ---   Stim animation interrupted   ---')
% %     else
% %         psychrethrow(psychlasterror);
% %     end
% % %     LoadIdentityClut(screen);
% % 	Screen('CloseAll');
% %     Priority(0);
% % 	ShowCursor;
% %
% %     % Bye bye!
% %     clearvars -global -except Ard
% %     clearvars;
% %     return;


% end % try-catch



end
%
function Drift_Simple_2
clearvars -global -except Ard; clearvars; sca;
% clear all; clear global; clear mex; sca; %#ok<CLMEX,CLALL>
[scriptPath, scriptName] = fileparts(mfilename('fullpath')); %#ok<*ASGLU>
TimeAnimation=0; stim_id_list=[]; 
Screen('Preference', 'SkipSyncTests', 1); %to skip synchronization failure set to 1

     calc_frames = 0; frame_time = 0.114;%0.113%0.0088%0.044%0.176%0.264;
     %define target eye (1 is left eye or 2 is right eye)
     tar_eye = 1;
     Param.tar_eye = tar_eye;
Param. Dichoptic = 1 ; Param.stereoMode = 0; %1
Param. screenid  = 0; % it matters only if Param.Dichoptic=0;
monitor_dist_cm  = 20;%13;
rot_bias_left   = -0 ; % if eye_rot_bias>0, set [+x,-x] L,R; else if eye_rot_bias<0, set [-x,+x] L,R
rot_bias_right  = +0 ;
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
    Param.StimProtocols.   DriftingGrating2          = 0 ; % set to 2 for "navigation stim"
    Param.StimProtocols.   DriftingGrating3          = 0 ;
	Param.StimProtocols.   DriftingGrating4          = 0 ;
    Param.StimProtocols.   DriftingGrating5          = 0 ;
	Param.StimProtocols.   DriftingGrating6          = 0 ;
    Param.StimProtocols.   DriftingGrating7          = 0 ;
    Param.StimProtocols.   DriftingGrating8          = 0 ;
    Param.StimProtocols.   DriftingGrating9          = 0 ;
    

    
%% DG
reps       = 4;
stim       = 5;
tf         = [1 2 4 1 2 4 1 2 4];
sf         = [0.02 0.02 0.02 0.04 0.04 0.04 0.08 0.08 0.08];
directions = 12;
postStim = 3;
angles     = [];%[0 90 180 270];
if Param.StimProtocols.DriftingGrating
    global DG %#ok<*TLEV>
	DG.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG.drawMask = 0; % draw grating through circular aperture
	DG.maskCenter_offset = [+0/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG.maskSize_deg = 70; %deg
	DG.spacfreq = sf(1); %cpd  
	DG.cyclespersecond = tf(1) ;
	DG.directions = directions; % Set the number of directions you want to be presented:
    DG.angles = angles;%[45 90 225 270];%[]; % specify cartesian angles
    DG.offset_rot_deg = [  rot_bias_left  ] ; %#ok<*NBRAK>
	DG.stimulus_time = stim ; % seconds per stimulus direction
	DG.poststim_time = postStim ; % seconds per blank stimulus (gray screen that precedes the stimulus)
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
    DG.time = DG.n_reps*DG.directions*(DG.poststim_time+DG.stimulus_time);
    TimeAnimation = TimeAnimation + DG.time;
end
if Param.StimProtocols.DriftingGrating2
    global DG2
	DG2.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG2.drawMask = 0; % draw grating through circular aperture
    DG2.maskCenter_offset = [+0/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG2.maskSize_deg = 70; %deg
	DG2.spacfreq = sf(2); %cpd  
	DG2.cyclespersecond = tf(2) ;
	DG2.directions = directions; % Set the number of directions you want to be presented:
    DG2.angles = angles;%[45 90 225 270];%[]; % specify cartesian angles
    DG2.offset_rot_deg = [ rot_bias_right ] ;
	DG2.stimulus_time = stim ; % seconds per stimulus direction
	DG2.poststim_time = postStim ;  % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG2.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG2.phase = [] ;
    DG2.amplitude = 0.35; % contrast (if 0.5->100%)
    DG2.BackgroundLuminance = 127;
    DG2.stim_id = 1.3; stim_id_list=[stim_id_list;DG2.stim_id];   % for the 'putsample' value
    if Param.StimProtocols.DriftingGrating2==2
        DG2.n_reps = 300 ;
        DG2.drawMask = 0;
        DG2.spacfreq = 0.04; %cpd  
        DG2.cyclespersecond = 2 ;
        DG2.directions = 12; DG2.angles = [];
        DG2.stimulus_time = 2 ;
        DG2.poststim_time = 3 ;
        Param.prestimInterval  = 1 ; %sec
    end
    DG2.time = DG2.n_reps*DG2.directions*(DG2.poststim_time+DG2.stimulus_time);
    TimeAnimation = TimeAnimation + DG2.time;
end
if Param.StimProtocols.DriftingGrating3
    global DG3
	DG3.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG3.drawMask = 0; % draw grating through circular aperture
    DG3.maskCenter_offset = [+3/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG3.maskSize_deg = 30; %deg
	DG3.spacfreq = sf(3); %cpd  
	DG3.cyclespersecond = tf(3) ;
	DG3.directions = directions; % Set the number of directions you want to be presented:
    DG3.angles = angles;%[]; % specify cartesian angles
	DG3.stimulus_time = stim ; % seconds per stimulus direction
	DG3.poststim_time = postStim ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG3.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG3.phase = [];
    DG3.amplitude = 0.35; % contrast (if 0.5->100%)
    DG3.offset_rot_deg = [  rot_bias_left ] ;
    DG3.BackgroundLuminance = 127;
    DG3.stim_id = 1.6;  stim_id_list=[stim_id_list;DG3.stim_id];  % for the 'putsample' value
    DG3.time = DG3.n_reps*DG3.directions*(DG3.poststim_time+DG3.stimulus_time);
    TimeAnimation = TimeAnimation + DG3.time;
end
if Param.StimProtocols.DriftingGrating4
    global DG4
	DG4.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG4.drawMask = 0; % draw grating through circular aperture
    DG4.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG4.maskSize_deg = 30; %deg
	DG4.spacfreq = sf(4); %cpd  
	DG4.cyclespersecond = tf(4) ;
	DG4.directions = directions; % Set the number of directions you want to be presented:
    DG4.angles = angles;%[]; % specify cartesian angles
	DG4.stimulus_time = stim ; % seconds per stimulus direction
	DG4.poststim_time = postStim ; % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG4.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG4.phase = [];
    DG4.amplitude = 0.35; % contrast (if 0.5->100%)
    DG4.offset_rot_deg = [  rot_bias_right ] ;
    DG4.BackgroundLuminance = 127;
    DG4.stim_id = 2;  stim_id_list=[stim_id_list;DG4.stim_id];  % for the 'putsample' value
    DG4.time = DG4.n_reps*DG4.directions*(DG4.poststim_time+DG4.stimulus_time);
    TimeAnimation = TimeAnimation + DG4.time;
end
if Param.StimProtocols.DriftingGrating5
    global DG5
	DG5.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG5.drawMask = 0; % draw grating through circular aperture
    DG5.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG5.maskSize_deg = 30; %deg
	DG5.spacfreq = sf(5); %cpd  
	DG5.cyclespersecond = tf(5) ;
	DG5.directions = directions; % Set the number of directions you want to be presented:
    DG5.angles = angles;%[]; % specify cartesian angles
	DG5.stimulus_time = stim ; % seconds per stimulus direction
	DG5.poststim_time = postStim ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG5.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG5.phase = [];
    DG5.amplitude = 0.35; % contrast (if 0.5->100%)
    DG5.offset_rot_deg = [  rot_bias_left  ] ;
    DG5.BackgroundLuminance = 127;
    DG5.stim_id = 2.3;  stim_id_list=[stim_id_list;DG5.stim_id];  % for the 'putsample' value
    DG5.time = DG5.n_reps*DG5.directions*(DG5.poststim_time+DG5.stimulus_time);
    TimeAnimation = TimeAnimation + DG5.time;
end
if Param.StimProtocols.DriftingGrating6
    global DG6
	DG6.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG6.drawMask = 0; % draw grating through circular aperture
    DG6.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG6.maskSize_deg = 30; %deg
	DG6.spacfreq = sf(6); %cpd  
	DG6.cyclespersecond = tf(6) ;
	DG6.directions = directions; % Set the number of directions you want to be presented:
    DG6.angles = angles;%[]; % specify cartesian angles
	DG6.stimulus_time = stim ; % seconds per stimulus direction
    DG6.poststim_time = postStim ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG6.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG6.phase = [];
    DG6.amplitude = 0.35; % contrast (if 0.5->100%)
    DG6.offset_rot_deg = [  rot_bias_right  ] ;
    DG6.BackgroundLuminance = 127;
    DG6.stim_id = 2.6;  stim_id_list=[stim_id_list;DG6.stim_id];  % for the 'putsample' value
    DG6.time = DG6.n_reps*DG6.directions*(DG6.poststim_time+DG6.stimulus_time);
    TimeAnimation = TimeAnimation + DG6.time;
end
if Param.StimProtocols.DriftingGrating7
    global DG7
	DG7.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG7.drawMask = 0; % draw grating through circular aperture
    DG7.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG7.maskSize_deg = 30; %deg
	DG7.spacfreq = sf(7); %cpd  
	DG7.cyclespersecond = tf(7) ;
	DG7.directions = directions; % Set the number of directions you want to be presented:
    DG7.angles = angles;%[]; % specify cartesian angles
	DG7.stimulus_time = stim ; % seconds per stimulus direction
    DG7.poststim_time = postStim ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG7.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG7.phase = [];
    DG7.amplitude = 0.35; % contrast (if 0.5->100%)
    DG7.offset_rot_deg = [  rot_bias_right  ] ;
    DG7.BackgroundLuminance = 127;
    DG7.stim_id = 3;  stim_id_list=[stim_id_list;DG7.stim_id];  % for the 'putsample' value
    DG7.time = DG7.n_reps*DG7.directions*(DG7.poststim_time+DG7.stimulus_time);
    TimeAnimation = TimeAnimation + DG7.time;
end
if Param.StimProtocols.DriftingGrating8
    global DG8
	DG8.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG8.drawMask = 0; % draw grating through circular aperture
    DG8.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG8.maskSize_deg = 30; %deg
	DG8.spacfreq = sf(8); %cpd  
	DG8.cyclespersecond = tf(8) ;
	DG8.directions = directions; % Set the number of directions you want to be presented:
    DG8.angles = angles;%[]; % specify cartesian angles
	DG8.stimulus_time = stim ; % seconds per stimulus direction
    DG8.poststim_time = postStim ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG8.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG8.phase = [];
    DG8.amplitude = 0.35; % contrast (if 0.5->100%)
    DG8.offset_rot_deg = [  rot_bias_right  ] ;
    DG8.BackgroundLuminance = 127;
    DG8.stim_id = 3.3;  stim_id_list=[stim_id_list;DG8.stim_id];  % for the 'putsample' value
    DG8.time = DG8.n_reps*DG8.directions*(DG8.poststim_time+DG8.stimulus_time);
    TimeAnimation = TimeAnimation + DG8.time;
end
if Param.StimProtocols.DriftingGrating9
    global DG9
	DG9.n_reps = reps ; % number of repetitions (1 repetition = each direction once).
	DG9.drawMask = 0; % draw grating through circular aperture
    DG9.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG9.maskSize_deg = 30; %deg
	DG9.spacfreq = sf(9); %cpd  
	DG9.cyclespersecond = tf(9) ;
	DG9.directions = directions; % Set the number of directions you want to be presented:
    DG9.angles = angles;%[]; % specify cartesian angles
	DG9.stimulus_time = stim ; % seconds per stimulus direction
    DG9.poststim_time = postStim ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG9.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG9.phase = [];
    DG9.amplitude = 0.35; % contrast (if 0.5->100%)
    DG9.offset_rot_deg = [  rot_bias_right  ] ;
    DG9.BackgroundLuminance = 127;
    DG9.stim_id = 3.6;  stim_id_list=[stim_id_list;DG9.stim_id];  % for the 'putsample' value
    DG9.time = DG9.n_reps*DG9.directions*(DG9.poststim_time+DG9.stimulus_time);
    TimeAnimation = TimeAnimation + DG9.time;
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
            Param.screenid(1) = 2 ;
            Param.screenid(2) = 3 ;%DRAGO

        end
            mouse_dist_cm([1 2]) = [monitor_dist_cm monitor_dist_cm] ;
            screenSize_cm_horz = 53;%44;
    elseif Param.stereoMode == 0
%         Param.screenid(1) = 1 ;%DRAGO
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
        mouse_dist_cm = 13%20%13.3%#ok<*NOPRT> %12.7%19 %18.8 ;
        screenSize_cm_horz = 44;%47.5%44;
    elseif Param.screenid == 2
        mouse_dist_cm = 12%20 %18.8 ;
        screenSize_cm_horz = 44;
    elseif Param.screenid == 0
        mouse_dist_cm = 20%18.8 %18.8 ;
        screenSize_cm_horz = 88;
    end
end

TimeAnimation = TimeAnimation + Param.prestimInterval;
    framesTOT = ceil( TimeAnimation / frame_time );
    disp(['Total time for whole stimulation is ' num2str(TimeAnimation/60) ' min.'])
    disp(['Total number of frames with frame time ' num2str(frame_time) ' is ' num2str(framesTOT) ])
    if calc_frames, return, end;
% Check if the stim_id of the stimuli are far enough:
if Param.useDAQdevice==1, thr_stimid=0.00001; else thr_stimid=0; end;
% if any(abs(diff(sort(stim_id_list))) <= thr_stimid )
%     warndlg({['The stim_id used are too close, you won''t be able to distinguish the different stimuli! Aborting...'];...
%              [ 'stim_id_list = ' mat2str(stim_id_list')]} )
%     return
% end

if     strfind(eval('computer'),'WIN')
    save_dir=['C:\Code\alessandro\stim_data\' datestr(now, 'yyyy_mm_dd')];
elseif strfind(eval('computer'),'LNX')
    save_dir=['/home/alessandro/Dropbox/Code/alessandro/stim_data/' datestr(now, 'yyyy_mm_dd')];
end
if ~isdir(save_dir), mkdir(save_dir), end;
datetime_suffix = datestr(now, 'yyyy_mm_dd_HH_MM_SS') ;
global stim_fname
stim_fname = [save_dir filesep 'InfoStim-' datetime_suffix '-DriftGrating_RetPos_PhiMotion'];
Param.save_dir = save_dir; Param.datetime_suffix = datetime_suffix;
Param.stim_fname = stim_fname;

internalRotation = 1; % 1=rotation of the gratings inside the "box", 0=rotation of the box.
if   internalRotation, Param.rotateMode = kPsychUseTextureMatrixForRotation;
else                   Param.rotateMode = []; end;
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
%             [screen(1), Param.screenRect ] = PsychImaging('OpenWindow',...
%             Param.screenid(1), 200, [],[],[], Param.stereoMode);%DRAGO
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
    pos_screen2 = [ mon_pos(1,1)+50  50  mon_pos(1,1)+50+round(mon_pos(1,3)*2/3)  round(mon_pos(1,4)*2/3) ];
    if Param.Dichoptic
        if Param.stereoMode == 4 || Param.stereoMode == 1
            [screen(1), Param.screenRect ] = PsychImaging('OpenWindow', Param.screenid(1), 200, pos_screen1,[],[], Param.stereoMode);
%             [screen(1), Param.screenRect ] = Screen('OpenWindow', Param.screenid(1), 127, pos_screen1,[],[], Param.stereoMode);
        elseif Param.stereoMode == 0
            [screen(1), Param.screenRect ] = Screen('OpenWindow', Param.screenid(1), 127, pos_screen1 );
            [screen(2), Param.screenRectR] = Screen('OpenWindow', Param.screenid(2), 127, pos_screen2 );
        end
    else %Param.Dichoptic==0
        [screen(1), Param.screenRect] = Screen('OpenWindow', Param.screenid(1), 127, pos_screen1 );
    end
end
if loadMyGammaTable == 0
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
if Param.Dichoptic
    degperpix = (2*atan(Param.pixel_cm./(2*mouse_dist_cm(2)))).*(180/pi);  %%true only for small angles, more generally would be: atan(pix./animal_dist_cm)  %*(180/pi) is necessary because atan in matlab gives result in radians
    Param.pixperdeg(2) = 1./degperpix;
end
Param.mouse_dist_cm = mouse_dist_cm;
Param.screenSize_cm_horz = screenSize_cm_horz;
%% Set parameters of stimulations
Param.stimSeq = {};
if Param.StimProtocols.DriftingGrating
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG, Param ] = SetDriftingGratings( screen(1), DG, Param )
end
if Param.StimProtocols.DriftingGrating2
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye);else win=screen(1); end
    if Param.StimProtocols.DriftingGrating
        DG2.BackgroundLuminance = DG.BackgroundLuminance ;
    end
    
    [ DG2, Param ] = SetDriftingGratings( win, DG2, Param )
end
if Param.StimProtocols.DriftingGrating3
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG3, Param ] = SetDriftingGratings( screen(1), DG3, Param )
end
if Param.StimProtocols.DriftingGrating4
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG4, Param ] = SetDriftingGratings( screen(1), DG4, Param )
end
if Param.StimProtocols.DriftingGrating5
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG5, Param ] = SetDriftingGratings( screen(1), DG5, Param )
end
if Param.StimProtocols.DriftingGrating6
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG6, Param ] = SetDriftingGratings( screen(1), DG6, Param )
end
if Param.StimProtocols.DriftingGrating7
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG7, Param ] = SetDriftingGratings( screen(1), DG7, Param )
end
if Param.StimProtocols.DriftingGrating8
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG8, Param ] = SetDriftingGratings( screen(1), DG8, Param )
end
if Param.StimProtocols.DriftingGrating9
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(tar_eye); else win=screen(1); end
    [ DG9, Param ] = SetDriftingGratings( screen(1), DG9, Param )
end

% make stimSeq a row numeric array:
Param.stimSeq = cat(2, Param.stimSeq{:} );
if Param.fullyInterleaved
	Param.stimSeq = Param.stimSeq(randperm(length(Param.stimSeq))); % elements 1 to 4 shuffled
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
       
if     Param.stereoMode==0
    for i=1:length(screen)
        Screen('FillRect', screen(i), prestimColor);
        Screen('Flip', screen(i));
    end
elseif Param.stereoMode==4 || Param.stereoMode == 1
    Screen('SelectStereoDrawBuffer', screen(1), 0);
    Screen('FillRect', screen(1), prestimColor);
    Screen('SelectStereoDrawBuffer', screen(1), 1);
    Screen('FillRect', screen(1), prestimColor);
    Screen('Flip', screen(1));
end
% Wait trigger to start stimulation protocol:
if test == 1 || Param.useDAQdevice==0 % Wait for release of all keys on keyboard, then sync us to retrace:
%     startKey = KbName('space'); keyCode = [];
% 	while ~keyCode(startKey)  %Ale: keep pausing as long as you don't press the startKey
%        [~,~,keyCode] = KbCheck; 
%        KbReleaseWait;
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
        [DG,  vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG ,tar_eye);
		    DG.cnt = DG.cnt+1;
    elseif exist('DG2','var') && Param.stimSeq(s)==DG2.stim_id % 
        [DG2, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG2 ,tar_eye);
        DG2.cnt = DG2.cnt+1;
    elseif exist('DG3','var') && Param.stimSeq(s)==DG3.stim_id % 
        [DG3, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG3 ,tar_eye);
        DG3.cnt = DG3.cnt+1;
    elseif exist('DG4','var') && Param.stimSeq(s)==DG4.stim_id % 
        [DG4, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG4 ,tar_eye);
        DG4.cnt = DG4.cnt+1;
    elseif exist('DG5','var') && Param.stimSeq(s)==DG5.stim_id % 
        [DG5, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG5 ,tar_eye);
        DG5.cnt = DG5.cnt+1;
    elseif exist('DG6','var') && Param.stimSeq(s)==DG6.stim_id % 
        [DG6, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG6 ,tar_eye);
        DG6.cnt = DG6.cnt+1;
    elseif exist('DG7','var') && Param.stimSeq(s)==DG7.stim_id % 
        [DG7, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG7 ,tar_eye);
        DG7.cnt = DG7.cnt+1;
    elseif exist('DG8','var') && Param.stimSeq(s)==DG8.stim_id %
        [DG8, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG8 ,tar_eye);
        DG8.cnt = DG8.cnt+1;
    elseif exist('DG9','var') && Param.stimSeq(s)==DG9.stim_id %
        [DG9, vbl] = PresentDriftingGrating_Simple( nidaq, screen, vbl, Param, DG9 ,tar_eye);
        DG9.cnt = DG9.cnt+1;
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
elseif Param.stereoMode==4 || Param.stereoMode == 1
    Screen('SelectStereoDrawBuffer', screen(1), 0);
    Screen('FillRect', screen(1), 80);
    Screen('SelectStereoDrawBuffer', screen(1), 1);
    Screen('FillRect', screen(1), 80);
    Screen('Flip', screen(1));
end

%% Closes windows and save parameters:

SaveSimple_2;
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

for screen_var = 1:length(screen)
    LoadIdentityClut(screen(screen_var));
end
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
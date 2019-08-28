%
function DriftGrating_RetPosition_PhiMotion
clearvars -global -except Ard; clearvars; sca;
% clear all; clear global; clear mex; sca; %#ok<CLMEX,CLALL>
[scriptPath, scriptName] = fileparts(mfilename('fullpath')); %#ok<*ASGLU>
TimeAnimation=0; stim_id_list=[]; 
Screen('Preference', 'SkipSyncTests', 1); %to skip synchronization failure set to 1

     calc_frames = 0 ; frame_time = 0.113;%0.113%0.0088%0.044%0.176%0.264;
Param. Dichoptic = 1 ; Param.stereoMode = 0; %1
Param. screenid  = 0; % it matters only if Param.Dichoptic=0;
monitor_dist_cm  = 13;%13;
rot_bias_left   = -0 ; % if eye_rot_bias>0, set [+x,-x] L,R; else if eye_rot_bias<0, set [-x,+x] L,R
rot_bias_right  = +0 ;
% Set random or sequential sequence of stimuli (directions or positions):
Param.seqmode='random'; %'sequential'
% Param.seqmode='sequential';
Param.fullyInterleaved = 1 ; %
Param.prestimInterval  = 0 ; %sec
Param.usePhotodiode = 1;
Param.useDAQdevice  = 0;
Param.useArduino    = 1;
test = 0 ;          % set to 1 if you want to test the stimulation without external triggering
test_screen = 0 ;   % set to 1 if you want to test the stimulation on a smaller screen; 0=full field.
loadMyGammaTable = 0 ;

    Param.StimProtocols.   DriftingGrating           = 1 ; % set to 2 for "navigation stim"
    Param.StimProtocols.   DriftingGrating2          = 1 ; % set to 2 for "navigation stim"
    Param.StimProtocols.   DriftingGrating3          = 0 ;
	Param.StimProtocols.   DriftingGrating4          = 0 ;
    Param.StimProtocols.   DriftingGrating5          = 0 ;
	Param.StimProtocols.   DriftingGrating6          = 0 ;
    Param.StimProtocols.DriftingGratingDisparity     = 0 ;
    Param.StimProtocols.DriftingGratingDisparity2    = 0 ;
    Param.StimProtocols.DriftingGratingDisparity3    = 0 ;
    Param.StimProtocols.DriftingGratingFlow          = 0 ;
    Param.StimProtocols.DriftingGratingFlow2         = 0 ;
    Param.StimProtocols.DriftingGratingFlow3         = 0 ;
    Param.StimProtocols. DriftingGratingInDepth         = 0 ;
    Param.StimProtocols.    DriftingGratingPhi       = 0 ;
    Param.StimProtocols.    StaticGrating            = 0 ;
    Param.StimProtocols. RetinotopicPosition         = 0 ;
	Param.StimProtocols.    RetinotopicPositionStripes  = 0 ;
    Param.StimProtocols. RFMapping                   = [0] ; %%[1,2]
    Param.StimProtocols.PhiMotionGratings            = 0 ;
    Param.StimProtocols.DriftingGratingDisparity_phi = 0 ;
    Param.StimProtocols.PhiMotionGratingsControl     = 0 ;
    Param.StimProtocols.  PhiGratings                = 0 ;
    Param.StimProtocols.     PhiMotionBars           = 0 ;
    Param.StimProtocols.        DotMotion           = 0 ;
    Param.StimProtocols.        DotMotion2          = 0 ;
    Param.StimProtocols.        DotMotion3          = 0 ;
    Param.StimProtocols.  RandomDotStereogram       = 0 ;
    Param.StimProtocols.  RandomDotStereogram2      = 0 ;
    Param.StimProtocols.  RandomBarStereogram       = 0 ;
    Param.StimProtocols.  RandomBarStereogram2      = 0 ;
    Param.StimProtocols.  NaturalImageStereogram    = 0 ;
    Param.StimProtocols.       RandomDotStereoBar   = 0 ;
    Param.StimProtocols.       RandomDotStereoEdge  = 0 ;

    
%% DG
reps       = [3 5 5];
stim       = [3 3 3];
tf         = [2 2 4];
sf         = [0.05 0.05 0.05];
directions = [16 4 4];
angles     = [];%[0 90 180 270];
if Param.StimProtocols.DriftingGrating
    global DG %#ok<*TLEV>
	DG.n_reps = reps(1) ; % number of repetitions (1 repetition = each direction once).
	DG.drawMask = 0; % draw grating through circular aperture
	DG.maskCenter_offset = [+0/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG.maskSize_deg = 70; %deg
	DG.spacfreq = sf(1); %cpd  
	DG.cyclespersecond = tf(1) ;
	DG.directions = directions(1); % Set the number of directions you want to be presented:
    DG.angles = angles;%[45 90 225 270];%[]; % specify cartesian angles
    DG.offset_rot_deg = [  rot_bias_left  ] ; %#ok<*NBRAK>
	DG.stimulus_time = stim(1) ; % seconds per stimulus direction
	DG.poststim_time = 1 ; % seconds per blank stimulus (gray screen that precedes the stimulus)
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
	DG2.n_reps = reps(1) ; % number of repetitions (1 repetition = each direction once).
	DG2.drawMask = 0; % draw grating through circular aperture
    DG2.maskCenter_offset = [+0/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG2.maskSize_deg = 70; %deg
	DG2.spacfreq = sf(1); %cpd  
	DG2.cyclespersecond = tf(1) ;
	DG2.directions = directions(1); % Set the number of directions you want to be presented:
    DG2.angles = angles;%[45 90 225 270];%[]; % specify cartesian angles
    DG2.offset_rot_deg = [ rot_bias_right ] ;
	DG2.stimulus_time = stim(1) ; % seconds per stimulus direction
	DG2.poststim_time = 1 ;  % seconds per blank stimulus (gray screen that precedes the stimulus)
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
	DG3.n_reps = reps(2) ; % number of repetitions (1 repetition = each direction once).
	DG3.drawMask = 0; % draw grating through circular aperture
    DG3.maskCenter_offset = [+3/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG3.maskSize_deg = 30; %deg
	DG3.spacfreq = sf(2); %cpd  
	DG3.cyclespersecond = tf(2) ;
	DG3.directions = directions(2); % Set the number of directions you want to be presented:
    DG3.angles = angles;%[]; % specify cartesian angles
	DG3.stimulus_time = stim(2) ; % seconds per stimulus direction
	DG3.poststim_time = 2 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
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
	DG4.n_reps = reps(2) ; % number of repetitions (1 repetition = each direction once).
	DG4.drawMask = 0; % draw grating through circular aperture
    DG4.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG4.maskSize_deg = 30; %deg
	DG4.spacfreq = sf(2); %cpd  
	DG4.cyclespersecond = tf(2) ;
	DG4.directions = directions(2); % Set the number of directions you want to be presented:
    DG4.angles = angles;%[]; % specify cartesian angles
	DG4.stimulus_time = stim(2) ; % seconds per stimulus direction
	DG4.poststim_time = 2 ; % seconds per blank stimulus (gray screen that precedes the stimulus)
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
	DG5.n_reps = reps(3) ; % number of repetitions (1 repetition = each direction once).
	DG5.drawMask = 0; % draw grating through circular aperture
    DG5.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG5.maskSize_deg = 30; %deg
	DG5.spacfreq = sf(3); %cpd  
	DG5.cyclespersecond = tf(3) ;
	DG5.directions = directions(3); % Set the number of directions you want to be presented:
    DG5.angles = angles;%[]; % specify cartesian angles
	DG5.stimulus_time = stim(3) ; % seconds per stimulus direction
	DG5.poststim_time = 2 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
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
	DG6.n_reps = reps(3) ; % number of repetitions (1 repetition = each direction once).
	DG6.drawMask = 0; % draw grating through circular aperture
    DG6.maskCenter_offset = [+2/10, -1/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DG6.maskSize_deg = 30; %deg
	DG6.spacfreq = sf(3); %cpd  
	DG6.cyclespersecond = tf(3) ;
	DG6.directions = directions(3); % Set the number of directions you want to be presented:
    DG6.angles = angles;%[]; % specify cartesian angles
	DG6.stimulus_time = stim(3) ; % seconds per stimulus direction
    DG6.poststim_time = 2 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DG6.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DG6.phase = [];
    DG6.amplitude = 0.35; % contrast (if 0.5->100%)
    DG6.offset_rot_deg = [  rot_bias_right  ] ;
    DG6.BackgroundLuminance = 127;
    DG6.stim_id = 2.6;  stim_id_list=[stim_id_list;DG6.stim_id];  % for the 'putsample' value
    DG6.time = DG6.n_reps*DG6.directions*(DG6.poststim_time+DG6.stimulus_time);
    TimeAnimation = TimeAnimation + DG6.time;
end
%% DGD
if Param.StimProtocols.DriftingGratingDisparity
    global DGD
	DGD.n_reps = 5 ; % number of repetitions (1 repetition = each direction once).
	DGD.drawMask = 0; % draw grating through circular aperture
	DGD.maskCenter_offset = [-3/10, 0/6 ; -3/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DGD.maskSize_deg = [80 ; 80]; %deg
    DGD.amplitude = 0.3;
	DGD.spacfreq = 0.01; %cpd  
	DGD.cyclespersecond = 2 ;
	DGD.directions = 2; % Set the number of directions you want to be presented:
    DGD.angles = [ 90 270 ];%[ 45 90 225 270];%[]; % specify cartesian angles
	DGD.stimulus_time = 2 ; % seconds per stimulus direction
    DGD.poststim_time = 2 ; % seconds per blank stimulus (gray screen that precedes the stimulus)
	DGD.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    % set interocular phase differences:
    DGD.nPhases = 12; % if 8, 0:45:315
    %set the center of RF (rightward, downward from screen center, values in % of screen):
    DGD.offset_phase1_perc = [ -0/18    -0/10 ];%[0 0];%[ 0.75/10 0.4/6 ];%[2/3 1/2];%[  2.5/10 2/3  ] ; 
    DGD.offset_phase2_perc = [ -0/18    -0/10 ];%[0 0];%[ 3.75/10 0.5/6 ] ;
    DGD.offset_rot_deg = [ rot_bias_left, rot_bias_right ] ;
    DGD.BackgroundLuminance = 127;
    Param.stereoMode = 1;
    DGD.stim_id = 3;  stim_id_list=[stim_id_list;DGD.stim_id];  % for the 'putsample' value
    DGD.time = DGD.n_reps*DGD.directions*DGD.nPhases*(DGD.poststim_time+DGD.stimulus_time);
    TimeAnimation = TimeAnimation + DGD.time;
end
if Param.StimProtocols.DriftingGratingDisparity2
    global DGD2
    DGD2 = DGD;
%     DGD2.n_reps = 3 ;
%     DGD2.offset_rot_deg = [ rot_bias_left, rot_bias_right] ;
	DGD2.spacfreq = 0.10; %cpd  
	DGD2.cyclespersecond = 2 ;
    DGD2.stim_id = 3.3;  stim_id_list=[stim_id_list;DGD2.stim_id];  % for the 'putsample' value
    DGD2.time = DGD2.n_reps*DGD2.directions*DGD2.nPhases*(DGD2.poststim_time+DGD2.stimulus_time);
    TimeAnimation = TimeAnimation + DGD2.time;
end
if Param.StimProtocols.DriftingGratingDisparity3
    global DGD3
    DGD3 = DGD;
%     DGD3.n_reps = 3;
	DGD3.spacfreq = 0.01; %cpd  
	DGD3.cyclespersecond = 2 ;
    DGD3.stim_id = 3.6;  stim_id_list=[stim_id_list;DGD3.stim_id];  % for the 'putsample' value
    DGD3.time = DGD3.n_reps*DGD3.directions*DGD3.nPhases*(DGD3.poststim_time+DGD3.stimulus_time);
    TimeAnimation = TimeAnimation + DGD3.time;
end
%% DGF
if Param.StimProtocols.DriftingGratingFlow
    global DGF
    if Param.StimProtocols.DriftingGratingDisparity
        DGF = DGD;
        DGF.n_reps = 2 ; % number of repetitions (1 repetition = each direction once).
        DGF.directions = 2; % Set the number of directions you want to be presented.
        DGF.angles = [90 270];%[]; % specify cartesian angles
    else
        DGF.n_reps = 4 ; % number of repetitions (1 repetition = each direction once).
        DGF.drawMask = 0; % draw grating through circular aperture
        DGF.maskCenter_offset = [+2/10, 1/6;
                                 +3/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
        DGF.maskSize_deg = [70;
                            70]; %deg
        DGF.spacfreq = 0.01; %cpd  
        DGF.cyclespersecond = 2 ;
        DGF.directions = 2; % Set the number of directions you want to be presented.
        DGF.angles = [ 90  270];%[]; % specify cartesian angles
        DGF.stimulus_time = 2 ; % seconds per stimulus direction
        DGF.poststim_time = 4 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
        DGF.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
        % set interocular phase differences:
        DGF.nPhases = 8; % if 8, 0:45:315
        DGF.offset_rot_deg = [ rot_bias_left, rot_bias_right] ;
        DGF.BackgroundLuminance = 127;
    end
    Param.stereoMode = 1;
    DGF.stim_id = 2;  stim_id_list=[stim_id_list;DGF.stim_id];  % for the 'putsample' value
    DGF.time = DGF.n_reps*DGF.directions*DGF.nPhases*(DGF.poststim_time+DGF.stimulus_time);
    TimeAnimation = TimeAnimation + DGF.time;
end
if Param.StimProtocols.DriftingGratingFlow2
    global DGF2
    DGF2 = DGF;
	DGF2.spacfreq = 0.01; %cpd  
	DGF2.cyclespersecond = 2 ;
    DGF2.stim_id = 4.3;  stim_id_list=[stim_id_list;DGF2.stim_id];  % for the 'putsample' value
    DGF2.time = DGF2.n_reps*DGF2.directions*DGF2.nPhases*(DGF2.poststim_time+DGF2.stimulus_time);
    TimeAnimation = TimeAnimation + DGF2.time;
end
if Param.StimProtocols.DriftingGratingFlow3
    global DGF3
    DGF3 = DGF;
	DGF3.spacfreq = 0.04; %cpd  
	DGF3.cyclespersecond = 2 ;
    DGF3.stim_id = 4.6;  stim_id_list=[stim_id_list;DGF3.stim_id];  % for the 'putsample' value
    DGF3.time = DGF3.n_reps*DGF3.directions*DGF3.nPhases*(DGF3.poststim_time+DGF3.stimulus_time);
    TimeAnimation = TimeAnimation + DGF3.time;
end
%% DGID
if Param.StimProtocols.DriftingGratingInDepth
    global DGID
	DGID.n_reps = 8 ; % number of repetitions (1 repetition = each direction once).
	DGID.drawMask = 0; % draw grating through circular aperture
	DGID.maskCenter_offset = [-3/10, 0/6 ; -3/10, 0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DGID.maskSize_deg = [80 ; 80]; %deg
    DGID.amplitude = 0.35;
	DGID.spacfreq = 0.01; %cpd  
	DGID.cyclespersecond = [1 2 4]; %Hz
	DGID.orientations = 1; % Set the number of directions you want to be presented:
    DGID.angles = [ 90 ];%[ 45 90 ];%[]; % specify cartesian angles
	DGID.stimulus_time = 2 ; % seconds per stimulus direction
    DGID.poststim_time = 1 ; % seconds per blank stimulus (gray screen that precedes the stimulus)
	DGID.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DGID.offset_rot_deg = [ rot_bias_left, rot_bias_right ] ;
    DGID.BackgroundLuminance = 127;
    DGID.stim_id = 4;  stim_id_list=[stim_id_list;DGID.stim_id];  % for the 'putsample' value
    DGID.time = DGID.n_reps*DGID.orientations*((numel(DGID.cyclespersecond)*2).^2)*(DGID.poststim_time+DGID.stimulus_time);
    TimeAnimation = TimeAnimation + DGID.time;
end
%% DGP & SG
if Param.StimProtocols.DriftingGratingPhi
    % Only two "frames" flashed.
    % DGP does not support dichoptic screens (the two "frames" are
    % presented to the same screen).
    global DGP
	DGP.n_reps = 4 ; % number of repetitions (1 repetition = each direction once).
	DGP.spacfreq = 0.04; %cpd  
	DGP.cyclespersecond = 0 ;
	DGP.directions = 8; % Set the number of directions you want to be presented:
    DGP.offset_rot_deg = [  0  ] ;
	DGP.poststim_time = 2.5 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DGP.patch_time = 0.125 ; % total seconds per stimulus direction (for each the "frames")
	DGP.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    DGP.phase = [0 45 90 135 180 225 270 315 ];%[ 0 45 90 135 ];[ 0 45 90 135 ];
    DGP.phaseincrement = 90 ;
    DGP.BackgroundLuminance = 127;
    DGP.stim_id = 2.3;  stim_id_list=[stim_id_list;DGP.stim_id];   % for the 'putsample' value
    DGP.time = DGP.n_reps*DGP.directions*(DGP.poststim_time+2*DGP.patch_time)*length(DGP.phase);
    TimeAnimation = TimeAnimation + DGP.time;
end
if Param.StimProtocols.StaticGrating
    % Static grating (same code as DGP but phaseincrement=0)
    % SG does not support dichoptic screens
    global SG
	SG.n_reps = 3 ; % number of repetitions (1 repetition = each orientation twice).
	SG.spacfreq = 0.04; %cpd  
	SG.cyclespersecond = 0 ;
	SG.directions = 8; % Set the number of directions you want to be presented:
    SG.offset_rot_deg = [  0  ] ;
	SG.poststim_time = 2.5 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	SG.stimulus_time = 0.125 ; % total seconds per stimulus direction
    SG.patch_time = SG.stimulus_time/2;
	SG.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    SG.phase = [0 45 90 135 180 225 270 315 ];%[ 0 45 90 135 ];
    SG.phaseincrement = 0 ;
    SG.BackgroundLuminance = 127;
    SG.stim_id = 2.6;  stim_id_list=[stim_id_list;SG.stim_id];   % for the 'putsample' value
    SG.time = SG.n_reps*SG.directions*(SG.poststim_time+2*SG.patch_time)*length(SG.phase);
    TimeAnimation = TimeAnimation + SG.time;
end
%% RP, RPS & RFM
if Param.StimProtocols.RetinotopicPosition
    global RP
	RP.n_reps = 3 ; % number of repetitions (1 repetition = each direction once).
	RP.spacfreq = 0.04; %cpd  
	RP.cyclespersecond = 2 ;
	RP.directions = 8; % Set the number of directions you want to be presented:
	RP.poststim_time = 2 ; % seconds per blank stimulus (gray screen that follows the stimulus)
	RP.stimulus_time = 4 ; % seconds per stimulus (all directions one after the other)
	RP.sinwave = 0;         % 1 = sine wave grating, 0 = square grating
    RP.BackgroundLuminance = 127;
	RP.stim_id = 0.2;  stim_id_list=[stim_id_list;RP.stim_id];
    RP.tot_positions = [5 3];%[10 4];%[10 8];%[5 4];%[10 8]; %[20 16];  % number of neighboring positions in X and Y where to present small grating stimuli in sequence.
    RP.margin.X = [0 0];%[0 1];%[4 2];%[6 6] ; % margin left and right
    RP.margin.Y = [0 0];%[0 0];%[2 2];%[7 0] ; % margin top and bottom
    RP.time = RP.n_reps*((RP.tot_positions(1)-sum(RP.margin.X))*(RP.tot_positions(2)-sum(RP.margin.Y)))*(RP.poststim_time+RP.stimulus_time);
    TimeAnimation = TimeAnimation + RP.time;
    if Param.Dichoptic && Param.stereoMode > 0;
        global RP2
        RP2=RP;
        RP2.n_reps = 2;%RP.n_reps ;
%         RP2.poststim_time = 0.1 ; % seconds per blank stimulus (gray screen that follows the stimulus)
%         RP2.stimulus_time = 0.1 ; % seconds per stimulus (all directions one after the other)
        RP2.margin.X = [0 0];%[4 2];%[6 6] ; % margin left and right
        RP2.margin.Y = [0 0];%[2 2];%[7 0] ; % margin top and bottom
        RP2.stim_id = 0.4;  stim_id_list=[stim_id_list;RP2.stim_id];
        RP2.time = RP2.n_reps*((RP2.tot_positions(1)-sum(RP2.margin.X))*(RP2.tot_positions(2)-sum(RP2.margin.Y)))*(RP2.poststim_time+RP2.stimulus_time);
        TimeAnimation = TimeAnimation + RP2.time;
    end
end
if Param.StimProtocols.RetinotopicPositionStripes
    global RPS
	RPS.n_reps = 8 ; %3 % number of repetitions (1 repetition = each direction once).
    RPS.amplitude = 0.5; %0.3
	RPS.spacfreq = 0.03; %cpd  
	RPS.cyclespersecond = 2 ;
	RPS.directions = 5; % Set the number of directions you want to be presented:
	RPS.poststim_time = 1.5 ; %1 % seconds per blank stimulus (gray screen that follows the stimulus)
	RPS.stimulus_time = 5 ; % seconds per stimulus (all directions one after the other)
	RPS.sinwave = 0;         % 1 = sine wave grating, 0 = square grating
    RPS.BackgroundLuminance = 127;
	RPS.stim_id = 0.6;  stim_id_list=[stim_id_list;RPS.stim_id];
    RPS.tot_positions = [11 6];%[10 4];%[10 8];%[5 4];%[10 8]; %[20 16];  % number of neighboring positions in X and Y where to present small grating stimuli in sequence.
    RPS.margin.X = [3 4];%[5 7];%[4 2];%[6 6] ; % margin left and right
    RPS.margin.Y = [0 0];%[0 0];%[2 2];%[7 0] ; % margin top and bottom
    RPS.OverlapStripes = 1; %if 1, stripes will overlap by 50%
    RPS.positions = [RPS.tot_positions(1)-sum(RPS.margin.X) , RPS.tot_positions(2)-sum(RPS.margin.Y)]; if RPS.OverlapStripes, RPS.positions = RPS.positions + (RPS.positions-1); end; RPS.positions(RPS.positions<0) = 0;
    RPS.time = RPS.n_reps*sum(RPS.positions)*(RPS.poststim_time+RPS.stimulus_time);
    TimeAnimation = TimeAnimation + RPS.time;
    if Param.Dichoptic && Param.stereoMode > 0
        global RPS2
        RPS2=RPS;
%         RPS2.n_reps =3 ;
%         RPS2.tot_positions = [18 6];
%         RPS2.margin.X = [6 6];%[4 2];%[6 6] ; % margin left and right
%         RPS2.margin.Y = [0 0];%[2 2];%[7 0] ; % margin top and bottom
        RPS2.stim_id = 0.8;  stim_id_list=[stim_id_list;RPS2.stim_id];
        RPS2.positions = [RPS2.tot_positions(1)-sum(RPS2.margin.X) , RPS2.tot_positions(2)-sum(RPS2.margin.Y)]; if RPS2.OverlapStripes, RPS2.positions = RPS2.positions + (RPS2.positions-1); end; RPS2.positions(RPS2.positions<0) = 0;
        RPS2.time = RPS2.n_reps*sum(RPS2.positions)*(RPS2.poststim_time+RPS2.stimulus_time);
        TimeAnimation = TimeAnimation + RPS2.time;
    end
end
if ismember(1,Param.StimProtocols.RFMapping)
    global RFM
    RFM.n_reps = 15 ;
    RFM.patch_time = 0.4;
    RFM.interpatch_time = 0.05;%0.030 ; % seconds bewteen patches.
    RFM.n_patches = [10, 9];
    RFM.fov = [80,90];  % to obtain patches of 10deg %[96,64];   %[96,80]; [84,72]  % size of the field of view in degree
    RFM.view_offset = [-17,-19];        % offset of the field of view (x<0 => sx , y<0 => up)
    RFM.rel_patch_size = 1.2;%1.2;           % patch size: 1: touching  - 0.5: size an distance is equal  %Ale: relative size, if 1 tha patches touch precisely, >1 they overlap
    RFM.Black = 1; RFM.BlackLuminance =  45;
    RFM.White = 1; RFM.WhiteLuminance = 210;
                   RFM.BackgroundLuminance = 127;
    RFM.gnomonic=0;
    RFM.stim_id = 3.6;  stim_id_list=[stim_id_list;RFM.stim_id];
    RFM.time = RFM.n_reps*prod(RFM.n_patches)*2*(RFM.patch_time+RFM.interpatch_time);
    TimeAnimation = TimeAnimation + RFM.time;
    if Param.StimProtocols.DriftingGrating || Param.StimProtocols.RetinotopicPosition || Param.StimProtocols.PhiGratings
        error('RFMapping must be the only StimProtocol ');
    end
    RFM.save_screen  = 1 ;
    show_screen_only = 0 ; % if 1 save and show the screen masks and exit without presenting any stimuli.
end
if ismember(2,Param.StimProtocols.RFMapping)
    global RFM2
    RFM2 = RFM;
%     RFM2.n_reps = 17 ;
    RFM2.stim_id = 3.8;  stim_id_list=[stim_id_list;RFM2.stim_id];
    RFM2.time = RFM2.n_reps*prod(RFM2.n_patches)*2*(RFM2.patch_time+RFM2.interpatch_time);
    TimeAnimation = TimeAnimation + RFM2.time;
    if Param.StimProtocols.DriftingGrating || Param.StimProtocols.RetinotopicPosition || Param.StimProtocols.PhiGratings
        error('RFMapping must be the only StimProtocol ');
    end
end
%% PMG
if Param.StimProtocols.PhiMotionGratings
    global PMG
	PMG.n_reps = 4 ; % number of repetitions (1 repetition = each direction once).
    PMG.amplitude = 0.25;
	PMG.spacfreq = 0.02; %cpd  
	PMG.cyclespersecond = 2;%2;%0.125 ; % TF of contrast reversal
	PMG.orientations = 2; % Set the number of orientations you want to be presented:
    PMG.angles = [0 90];%[ 45 90];%[]; % specify cartesian angles
    PMG.offset_rot_deg = [ rot_bias_left, rot_bias_right ] ; % angle offset of screen2 compared to screen1 (clockwise)
    %set the center of RF (rightward, downward from screen center, values in % of screen):
    PMG.offset_phase1_perc = [ -0/18    -0/10 ];%[0 0];%[ 0.75/10 0.4/6 ];%[2/3 1/2];%[  2.5/10 2/3  ] ; 
    PMG.offset_phase2_perc = [ -0/18    -0/10 ];%[0 0];%[ 3.75/10 0.5/6 ] ;
    PMG.stimulus_time = 3 ;
%     PMG.interpatch_time = 0;
	PMG.poststim_time = 1 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
% 	PMG.patch_time = 0.125 ; % total seconds per stimulus direction (for both the "frames")
	PMG.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    PMG.phases = [0 30 60 90 120 150 180 210 240 270 300 330];%[0 45 90 135 180 225 270 315];
%     PMG.phaseincrement = 90 ;
    PMG.BackgroundLuminance = 127;
    PMG.stim_id = 4;  stim_id_list=[stim_id_list;PMG.stim_id];   % for the 'putsample' value
    PMG.time = PMG.n_reps*PMG.orientations*(PMG.poststim_time+PMG.stimulus_time)*(length(PMG.phases)/2*length(PMG.phases));
    TimeAnimation = TimeAnimation + PMG.time;
end
%% DGDphi
if Param.StimProtocols.DriftingGratingDisparity_phi
    global DGDphi
	DGDphi.n_reps = 4 ; % number of repetitions (1 repetition = each direction once).
	DGDphi.drawMask = 0; % draw grating through circular aperture
	DGDphi.maskCenter_offset = [+2/10,1/6 ; +3/10,0/6]; %set the aperture center (rightward, downward from screen center, values in % of screen):
	DGDphi.maskSize_deg = [70 ; 70]; %deg
	DGDphi.spacfreq = 0.03; %cpd  
	DGDphi.cyclespersecond = 2 ;
	DGDphi.directions = 2; % Set the number of directions you want to be presented:
    DGDphi.angles = [90 270];%[0 90 180 270];%[ 45 90 225 270];%[]; % specify cartesian angles
	DGDphi.stimulus_time = 2.5 ; % seconds per stimulus direction
    DGDphi.poststim_time = 3 ; % seconds per blank stimulus (gray screen that precedes the stimulus)
	DGDphi.sinwave = 1;        % 1 = sine wave grating, 0 = square grating
    %set the center of RF (rightward, downward from screen center, values in % of screen):
    DGDphi.offset_phase1_perc = [ +1/8 -0.5/5 ];%[0 0];%[ 0.75/10 0.4/6 ];%[2/3 1/2];%[  2.5/10 2/3  ] ; 
    DGDphi.offset_phase2_perc = [ +1.25/8 -0.5/5 ];%[0 0];%[ 3.75/10 0.5/6 ] ;
    % set interocular phase differences:
%     DGDphi.nPhases = 8; % if 8, 0:45:315
    DGDphi.phases = [0:45:315];
    DGDphi.nPhases = numel(DGDphi.phases);
    DGDphi.offset_rot_deg = [ rot_bias_left, rot_bias_right ] ;
    DGDphi.BackgroundLuminance = 127;
    Param.stereoMode = 1;
    DGDphi.stim_id = 3.5;  stim_id_list=[stim_id_list;DGDphi.stim_id];  % for the 'putsample' value
    DGDphi.time = DGDphi.n_reps*DGDphi.directions*DGDphi.nPhases*2*(DGDphi.poststim_time+DGDphi.stimulus_time);
    TimeAnimation = TimeAnimation + DGDphi.time;
end

if Param.StimProtocols.PhiMotionGratingsControl
    global PMGc1 PMGc2
    % if Dichoptic==1 first (left) screen switches from phase 0 to blank to
    % 180 deg and second (right) screen switches from phase 90 to blank to
    % 270 deg, alternatively to produce phi motion "binocularly".
    % if Dichoptic==0 only first (left) screen performs phi motion,
    % switching from 0 to 90 to 180 to 270... without blank in between.
	PMGc1.n_reps = 1;%PMG.n_reps ; % number of repetitions (1 repetition = each direction once).
	PMGc1.spacfreq = PMG.spacfreq; %cpd  
	PMGc1.cyclespersecond = PMG.cyclespersecond;%2;%0.125 ;
	PMGc1.directions = PMG.directions; % Set the number of directions you want to be presented:
    PMGc1.offset_rot_deg = PMG.offset_rot_deg ; % angle offset of screen2 compared to screen1 (clockwise)
    %set the center of RF (rightward, downward from screen center, values in % of screen):
    PMGc1.offset_phase1_perc = PMG.offset_phase1_perc;%[2/3 1/2];%[  2.5/10 2/3  ] ; 
%     PMGc1.offset_phase2_perc = PMG.offset_phase2_perc ;
    PMGc1.stimulus_time = PMG.stimulus_time ;
    PMGc1.interpatch_time = PMG.interpatch_time;
	PMGc1.poststim_time = PMG.poststim_time ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
% 	PMG.patch_time = 0.125 ; % total seconds per stimulus direction (for both the "frames")
	PMGc1.sinwave = PMG.sinwave;        % 1 = sine wave grating, 0 = square grating
    PMGc1.phase = PMG.phase; %[180 225 270 315]; %[ 0 45 90 135];
    PMGc1.phaseincrement = PMG.phaseincrement ;
    PMGc1.BackgroundLuminance = PMG.BackgroundLuminance;
    PMGc1.stim_id = 4.2;  stim_id_list=[stim_id_list;PMGc1.stim_id];  % for the 'putsample' value
    PMGc1.time = PMGc1.n_reps*PMGc1.directions*(PMGc1.poststim_time+PMGc1.stimulus_time)*length(PMGc1.phase);
    %
    PMGc2 = PMGc1;
    PMGc2.offset_phase2_perc = PMG.offset_phase2_perc ;
    PMGc2 = rmfield(PMGc2, 'offset_phase1_perc') ;
    PMGc2.stim_id = 4.4;  stim_id_list=[stim_id_list;PMGc2.stim_id];
    TimeAnimation = TimeAnimation + 2* PMGc1.time;
end
%% PG
if Param.StimProtocols.PhiGratings
    global PG
    PG.n_reps = 1 ;
	PG.spacfreq = 0.04; %cpd  
    PG.sinwave = 1;
% 	PG.cyclespersecond = 2;%2;%0.125;  % PM.patch_time = 1/PM.cyclespersecond/4; % (0.125s for TF=2Hz)
    PG.phases = [ 0    72   144   216   288 ]; %[180 225 270 315]; %[ 0 45 90 135];
    PG.poststim_time = 2.5 ; % seconds per blank stimulus (gray screen that follows the stimulus)
	PG.patch_time = 0.250 ; % seconds per single flashed bar
    % seconds bewteen the first flashed grating and the second one (it's
    % not possible to set 0, the minimum will be one refresh rate 0.017):
    PG.interpatch_time = 0.017; 
    PG.BackgroundLuminance = 127;
    PG.directions = 8 ;
    PG.offset_rot_deg = [ 0  -33 ] ; % angle offset of screen1 and screen2 (clockwise)
    PG.stim_id = 4;  stim_id_list=[stim_id_list;PG.stim_id];
    PG.nbr_phases = length(PG.phases);
    PG.nbr_stimuli = (PG.nbr_phases^2)*(PG.directions/2)*2*2 ;  %*2 for dichoptic and monocular; *2 for screensequence.
    PG.time = PG.n_reps*PG.nbr_stimuli*(PG.patch_time*2+PG.interpatch_time+PG.poststim_time);
    TimeAnimation = TimeAnimation + PG.time;
end
%% PMB
if Param.StimProtocols.PhiMotionBars
    global PMB
    PMB.n_reps = 5 ;
    PMB.poststim_time = 1 ; % seconds per blank stimulus (gray screen that follows the stimulus)
	PMB.patch_time = 0.30 ; % seconds per single flashed bar
    PMB.interpatch_time = 0.017 ; % seconds bewteen the first flashed bar and the second one.
    PMB.Black = 0; PMB.BlackLuminance =   0;
    PMB.White = 1; PMB.WhiteLuminance = 255 - 30;
                  PMB.BackgroundLuminance = 127;
    PMB.fov = [140,120] ;
    % set offset of the bars from the centre of the FOV
    % (positive values -> right,bottom)
    % roughly: x: 1/10 screen=12.5, y: 1/6=11
    PMB.offset_stim(1,:) = [  -14  -10  ] ; 
    if Param.Dichoptic
%         PMB.offset_stim(2,:) = [  25.5  0  ] ;
        PMB.offset_stim(2,:) = PMB.offset_stim(1,:);
    end
    PMB.widthBar_deg = 10 ;
    PMB.lengthBar_deg = 120 ; %deg
    PMB.edgeSmoothing = 0 ;
    PMB.gnomonic = 0 ; % 1 if you want to use gnomonic projections, 0 not.
%     PMB.rotation_deg = [  0    45   90   135  ] ;
    PMB.rotation_deg = [     0     ] ;
    PMB.offset_rot_deg = [ rot_bias_left  rot_bias_right ] ;
    PMB.bar_centers = [0 5 10 15 20]-20/2;%[ 0 7 14 21 ] - 21/2;%[   -15 -5 +5 15    ] ;
%     PMB.bar_centers = [   -8 0 8    ] ;
    PMB.stim_id = 4;  stim_id_list=[stim_id_list;PMB.stim_id];
    if Param.Dichoptic
        PMB.stim_id = 4.5;
    end
    PMB.locations = length(PMB.bar_centers) ;
    PMB.nbr_stimuli = (PMB.Black+PMB.White)*(PMB.locations^2)*numel(PMB.rotation_deg)*2*2 ; %*2 for dichoptic and monocular; *2 for screensequence.
%     if Param.Dichoptic
%         PMB.nbr_stimuli = PMB.nbr_stimuli *2;
%     end
    PMB.time = PMB.n_reps*PMB.nbr_stimuli*(PMB.patch_time*2+PMB.interpatch_time+PMB.poststim_time);
    TimeAnimation = TimeAnimation + PMB.time;
    show_screen_only = 0; % if 1 save and show the screen masks and exit without presenting any stimuli.
    PMB.save_screen = 1;
%     if show_screen_only == 1; PMB.save_screen = 0;
%     else calc_frames = 0; end;
end
%% DM DotMotion
if Param.StimProtocols.DotMotion
    global DM
    DM.n_reps  = 6;
    DM.nDots   = []; % number of dots; if empty nDots is calculated from the density.
    DM.density = 0.20; % overall density, used only when nDots=[].
    DM.color{1} = [255,255,255]; % color of the dots
%     DM.color{2} = [0,0,0]; % if uncommented: hald dots are color{1} and half are color{2}.
    DM.BackgroundLuminance = 127;
    % diameter of dots (deg); with mouse_dist_cm=20 max size with 800x600 is 10, with 1280x720 is
    % 6, 
    DM.size = 4; 
    DM.apertureSize = []; % size of rectangular aperture [w,h] in degrees; if [], aperture=screen.
    DM.speed = 50;       %deg/sec (2Hz grating -> 50deg/s)
    DM.offset_rot_deg = [  0  ] ; %#ok<*NBRAK>
	DM.poststim_time = 3 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
	DM.stimulus_time = 3 ;
    DM.directions = 8;  %degrees (clockwise from straight up)
    DM.lifetime = 5;  %lifetime of each dot (sec)
    DM.stim_id = 2;  stim_id_list=[stim_id_list;DM.stim_id];   % for the 'putsample' value
    DM.time = DM.n_reps*DM.directions*(DM.poststim_time+DM.stimulus_time);
    TimeAnimation = TimeAnimation + DM.time;
end
if Param.StimProtocols.DotMotion2
    global DM2
    DM2=DM;
    DM2.nDots = []; % number of dots
    DM2.density = 0.10; % overall density, used only when nDots=[].
    DM2.size = 4; % size of dots (deg)
    DM2.speed = 50;       %deg/sec (2Hz grating -> 50deg/s)
    DM2.lifetime = DM2.stimulus_time;  %lifetime of each dot (sec)
    DM2.stim_id = 2.6; stim_id_list=[stim_id_list;DM2.stim_id];  % for the 'putsample' value
    DM2.time = DM2.n_reps*DM2.directions*(DM2.poststim_time+DM2.stimulus_time);
    TimeAnimation = TimeAnimation + DM2.time;
end
if Param.StimProtocols.DotMotion3
    global DM3
    DM3=DM;
    DM3.nDots = []; % number of dots
    DM3.density = 0.20; % overall density, used only when nDots=[].
    DM3.size = 2; % size of dots (deg)
    DM3.speed = 50;       %deg/sec (2Hz grating -> 50deg/s)
    DM3.lifetime = DM3.stimulus_time;  %lifetime of each dot (sec)
    DM3.stim_id = 0.5;  stim_id_list=[stim_id_list;DM3.stim_id]; % for the 'putsample' value
    DM3.time = DM3.n_reps*DM3.directions*(DM3.poststim_time+DM3.stimulus_time);
    TimeAnimation = TimeAnimation + DM3.time;
end
%% RDS
if Param.StimProtocols.RandomDotStereogram
    global RDS
    RDS.n_reps = 9;
    RDS.nDots   = []; % number of dots; if empty nDots is calculated from the density.
    RDS.density = 0.25; % overall density, used only when nDots=[].
    RDS.size = 12; % size of dots (deg)
    RDS.disparity_range = [-30 30]; %deg
    RDS.n_disparity_intervals = 23; % better odd number
    RDS.pattern_time = 0.150 ; %sec
	RDS.stimulus_time = 5 ;
	RDS.poststim_time = 2 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
    RDS.orientations = 1;  %degrees (clockwise from straight up)
    RDS.angles = [ 90 ]; % cartesian angles
    RDS.color{1} = [255,255,255]-30; % color of the dots
%     RDS.color{2} = [0 0 0]+30; %[0,0,0]; % if not empty: hald dots are color{1} and half are color{2}.
    RDS.Anticorrelated = 0;
    RDS.BackgroundLuminance = 127;
    RDS.apertureSize = []; % size of rectangular aperture [w,h] in degrees; if [], aperture=screen.
%     RDS.offset_rot_deg = [  0  ] ; %not supported yet! %#ok<*NBRAK>
    RDS.lifetime = RDS.stimulus_time;%not supported! %lifetime of each dot (sec) 
    RDS.stim_id = 2;  stim_id_list=[stim_id_list;RDS.stim_id];   % for the 'putsample' value
    RDS.time = RDS.n_reps*RDS.orientations*RDS.n_disparity_intervals*(RDS.poststim_time+RDS.stimulus_time);
    TimeAnimation = TimeAnimation + RDS.time;
end
if Param.StimProtocols.RandomDotStereogram2
    global RDS2
    RDS2.n_reps = 6;
    RDS2.nDots   = []; % number of dots; if empty nDots is calculated from the density.
    RDS2.density = 0.25; % overall density, used only when nDots=[].
    RDS2.size = 12; % size of dots (deg)
    RDS2.disparity_range = [-30 30]; %deg
    RDS2.n_disparity_intervals = 19; % better odd number
    RDS2.pattern_time = 0.125 ; %sec
	RDS2.stimulus_time = 5 ;
	RDS2.poststim_time = 2 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
    RDS2.orientations = 1;  %degrees (clockwise from straight up)
    RDS2.angles = [90 ]; % cartesian angles
    RDS2.color{1} = [255,255,255]-30; % color of the dots
%     RDS2.color{2} = [0 0 0]+30; %[0,0,0]; % if not empty: hald dots are color{1} and half are color{2}.
    RDS2.Anticorrelated = 1;
    RDS2.BackgroundLuminance = 127;
    RDS2.apertureSize = []; % size of rectangular aperture [w,h] in degrees; if [], aperture=screen.
%     RDS2.offset_rot_deg = [  0  ] ; %not supported yet! %#ok<*NBRAK>
    RDS2.lifetime = RDS2.stimulus_time;%not supported! %lifetime of each dot (sec) 
    RDS2.stim_id = 2.5;  stim_id_list=[stim_id_list;RDS2.stim_id];   % for the 'putsample' value
    RDS2.time = RDS2.n_reps*RDS2.orientations*RDS2.n_disparity_intervals*(RDS2.poststim_time+RDS2.stimulus_time);
    TimeAnimation = TimeAnimation + RDS2.time;
end
%% RBS
if Param.StimProtocols.RandomBarStereogram
    global RBS
    RBS.n_reps = 6;
    RBS.nBars = [3];
    RBS.barWidth = [6 12];
    RBS.disparity_range = [-30 30]; %deg
    RBS.n_disparity_intervals = 15; % better odd number
    RBS.color{1} = [255,255,255]-30; % color of the dots
%     RBS.color{2} = [0,0,0]; % if not empty: hald dots are color{1} and half are color{2}.
    RBS.BackgroundLuminance = 127;
    RBS.apertureSize = []; % size of rectangular aperture [w,h] in degrees; if [], aperture=screen.
    RBS.pattern_time = 0.160 ; %sec
	RBS.stimulus_time = 5 ;
    RBS.poststim_time = 3 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
    RBS.orientations = 1;
    RBS.angles = [90 ]; % cartesian angles
    RBS.offset_rot_deg = [ rot_bias_left, rot_bias_right ] ;
    RBS.lifetime = RBS.stimulus_time;  %lifetime of each dot (sec)
    RBS.stim_id = 3.8;  stim_id_list=[stim_id_list;RBS.stim_id];   % for the 'putsample' value
    RBS.time = RBS.n_reps*RBS.orientations*RBS.n_disparity_intervals*(RBS.poststim_time+RBS.stimulus_time);
    TimeAnimation = TimeAnimation + RBS.time;
end
if Param.StimProtocols.RandomBarStereogram2
    global RBS2
    RBS2.n_reps = 6;
    RBS2.nBars = [3];
    RBS2.barWidth = [6 12];
    RBS2.disparity_range = [-30 30]; %deg
    RBS2.n_disparity_intervals = 15; % better odd number
    RBS2.color{1} = [255,255,255]-30; % color of the dots
%     RBS2.color{2} = [0,0,0]; % if not empty: hald dots are color{1} and half are color{2}.
    RBS2.BackgroundLuminance = 127;
    RBS2.apertureSize = []; % size of rectangular aperture [w,h] in degrees; if [], aperture=screen.
    RBS2.pattern_time = 0.080 ; %sec
	RBS2.stimulus_time = 5 ;
    RBS2.poststim_time = 3 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
    RBS2.orientations = 1;
    RBS2.angles = [90 ]; % cartesian angles
    RBS2.offset_rot_deg = [ rot_bias_left, rot_bias_right ] ;
    RBS2.lifetime = RBS2.stimulus_time;  %lifetime of each dot (sec)
    RBS2.stim_id = 4.2;  stim_id_list=[stim_id_list;RBS2.stim_id];   % for the 'putsample' value
    RBS2.time = RBS2.n_reps*RBS2.orientations*RBS2.n_disparity_intervals*(RBS2.poststim_time+RBS2.stimulus_time);
    TimeAnimation = TimeAnimation + RBS2.time;
end
%% NIS
if Param.StimProtocols.NaturalImageStereogram
    global NIS
    NIS.n_reps = 6;
    NIS.disparity_range = [-30 30]; %deg
    NIS.n_disparity_intervals = 19; % better odd number
    NIS.pattern_time = 0.100 ; %sec
	NIS.stimulus_time = 5 ;
	NIS.poststim_time = 2 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
    NIS.orientations = 1;  %degrees (clockwise from straight up)
    NIS.angles = [ 90 ]; % cartesian angles
    NIS.BackgroundLuminance = 127;
    NIS.stim_id = 4.5;  stim_id_list=[stim_id_list;NIS.stim_id];   % for the 'putsample' value
    NIS.time = NIS.n_reps*NIS.orientations*NIS.n_disparity_intervals*(NIS.poststim_time+NIS.stimulus_time);
    TimeAnimation = TimeAnimation + NIS.time;
end
%% RDSB
if Param.StimProtocols. RandomDotStereoBar
    global RDSB
    % Check computation and drawing of dots with tictoc: it must be less
    % than pattern_time. If not, reduce screen resolution or dot density,
    % or increase dot size:
    RDSB.test_RDSB_ComputationTime = 0;
    RDSB.test_RDSB_pattern = 0;  if RDSB.test_RDSB_pattern==1, warning('RDBS is running in test mode'); pause(2); end;
    RDSB.n_reps = 5;
    RDSB.n_cycles = 4;
%     RDSB.disparity_range = [-30 30]; %deg
    RDSB.disparity_values    = [ -4 0 +4   ]; %neg->near
    RDSB.density = 0.5; % overall density, used only when nDots=[].
    RDSB.size = 1; % size of dots (deg)
    RDSB.color{1} = 150; % color of the dots
%     RDS.color{2} = [0 0 0]+30; %[0,0,0]; % if not empty: hald dots are color{1} and half are color{2}.
    RDSB.Anticorrelated = 0;
    RDSB.BackgroundLuminance = 70;
%   RDSB.spacfreq = 0.05; %cpd  
% 	RDSB.cyclespersecond = 2 ;
% 	RDSB.stimulus_time   = 5 ;
    RDSB.pattern_time    = 0.050 ; %sec (refresh rate=8.3ms
	RDSB.poststim_time   = 1 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
    RDSB.directions      = 4; % Set the number of directions you want to be presented:
    RDSB.angles          = [  ];%[ 45 90 225 270];%[]; % specify cartesian angles
    RDSB.fov             = [ 130, 130];
    RDSB.view_offset     = [ -15, 0];
    RDSB.bar_width  =    16; %deg
    RDSB.bar_length =    140; %deg
    RDSB.bar_speed  =    40; %deg/s (40 deg/s==2Hz 0.05 cpd)
%     RDSB.pixel_density = 0.5; % fraction bright pixels
%     RDSB.white = 255-30; % set brightness of white
%     RDSB.black =   0+30; % set brightness of black
%     RDSB.BackgroundLuminance = round( (RDSB.white*RDSB.pixel_density+RDSB.black*(1-RDSB.pixel_density)) );
    RDSB.stim_id = 2;  stim_id_list=[stim_id_list;RDSB.stim_id];   % for the 'putsample' value
    RDSB.n_disparities = numel(RDSB.disparity_values);
    RDSB.stimulus_time_perCycle = RDSB.fov(1)/RDSB.bar_speed;
    RDSB.stimulus_time = RDSB.stimulus_time_perCycle *RDSB.n_cycles;
    RDSB.time = RDSB.n_reps*RDSB.directions*RDSB.n_disparities*(RDSB.poststim_time+RDSB.stimulus_time);
    TimeAnimation = TimeAnimation + RDSB.time;
end
%% RDSE
if Param.StimProtocols. RandomDotStereoEdge
    global RDSE
    % Check computation and drawing of dots with tictoc: it must be less
    % than pattern_time. If not, reduce screen resolution or dot density,
    % or increase dot size:
    RDSE.test_RDSE_ComputationTime = 0;
    RDSE.n_reps   = 5;
    RDSE.n_cycles = 4;
%     RDSE.disparity_range = [-30 30]; %deg
%     RDSE.n_disparity_intervals = 19; % better odd number
    RDSE.disparity_values    = [ 0 4 8 ]; % positive values create gap at right screen
    RDSE.density = 0.9; % overall density, used only when nDots=[].
    RDSE.size = 8; % size of dots (deg)
    RDSE.color{1} = 150; % color of the dots
%     RDS.color{2} = [0 0 0]+30; %[0,0,0]; % if not empty: hald dots are color{1} and half are color{2}.
    RDSE.Anticorrelated = 0;
    RDSE.BackgroundLuminance = 70;
% 	RDSE.stimulus_time   = 5 ;
    RDSE.pattern_time    = 0.050 ; %sec (refresh rate=8.3ms
	RDSE.poststim_time   = 0 ;    % seconds per blank stimulus (gray screen that precedes the stimulus)
    RDSE.directions      = 4; % Set the number of directions you want to be presented:
    RDSE.angles          = [ ];%[ 45 90 225 270];%[]; % specify cartesian angles
    RDSE.fov             = [ 130, 130]; %deg
    RDSE.view_offset     = [ -15, 0];
    RDSE.bar_speed  =    40; %deg/s (40 deg/s==2Hz 0.05 cpd)
%     RDSE.pixel_density = 0.5; % fraction bright pixels
%     RDSE.white = 255-30; % set brightness of white
%     RDSE.black =   0+30; % set brightness of black
%     RDSE.BackgroundLuminance = round( (RDSE.white*RDSE.pixel_density+RDSE.black*(1-RDSE.pixel_density)) );
    RDSE.stim_id = 2.5;  stim_id_list=[stim_id_list;RDSE.stim_id];   % for the 'putsample' value
    RDSE.n_disparities = numel(RDSE.disparity_values);
    RDSE.stimulus_time_perCycle = RDSE.fov(1)/RDSE.bar_speed;
    RDSE.stimulus_time = RDSE.stimulus_time_perCycle *RDSE.n_cycles;
    RDSE.time = RDSE.n_reps*RDSE.directions*RDSE.n_disparities*(RDSE.poststim_time+RDSE.stimulus_time);
    TimeAnimation = TimeAnimation + RDSE.time;
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
        screenSize_cm_horz = 53;%44;
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
if Param.useDAQdevice==1, thr_stimid=0.29; else thr_stimid=0; end;
if any(abs(diff(sort(stim_id_list))) <= thr_stimid )
    warndlg({['The stim_id used are too close, you won''t be able to distinguish the different stimuli! Aborting...'];...
             [ 'stim_id_list = ' mat2str(stim_id_list')]} )
    return
end

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
if test_screen==0 && ( Param.StimProtocols.RandomDotStereoEdge || Param.StimProtocols.RandomDotStereoBar )
    oldRes = Screen('Resolution', Param.screenid(1));
    if oldRes.width ~= 1600 || oldRes.height ~= 900
        disp(' '); disp(' ');
        disp('  Warning:  Changing screen resolution to 1600x900'); pause(2);
        oldRes = SetResolution(Param.screenid(1), 1600, 900);
    end
end
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
    if Param.StimProtocols.RetinotopicPosition
        DG.BackgroundLuminance = RP.BackgroundLuminance ;
    end
    [ DG, Param ] = SetDriftingGratings( screen(1), DG, Param )
end
if Param.StimProtocols.DriftingGrating2
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(2);else win=screen(1); end
    if Param.StimProtocols.DriftingGrating
        DG2.BackgroundLuminance = DG.BackgroundLuminance ;
    end
    if Param.StimProtocols.RetinotopicPosition
        DG2.BackgroundLuminance = RP.BackgroundLuminance ;
    end
    [ DG2, Param ] = SetDriftingGratings( win, DG2, Param )
end
if Param.StimProtocols.DriftingGrating3
    [ DG3, Param ] = SetDriftingGratings( screen(1), DG3, Param )
end
if Param.StimProtocols.DriftingGrating4
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(2); else win=screen(1); end
    [ DG4, Param ] = SetDriftingGratings( screen(1), DG4, Param )
end
if Param.StimProtocols.DriftingGrating5
    [ DG5, Param ] = SetDriftingGratings( screen(1), DG5, Param )
end
if Param.StimProtocols.DriftingGrating6
    if Param.Dichoptic && Param.stereoMode==0
        win=screen(2); else win=screen(1); end
    [ DG6, Param ] = SetDriftingGratings( screen(1), DG6, Param )
end
if Param.StimProtocols.DriftingGratingDisparity
    [ DGD, Param ] = SetDriftingGratingsDisparity( screen, DGD, Param )
end
if Param.StimProtocols.DriftingGratingDisparity2
    [ DGD2, Param ] = SetDriftingGratingsDisparity( screen, DGD2, Param )
end
if Param.StimProtocols.DriftingGratingDisparity3
    [ DGD3, Param ] = SetDriftingGratingsDisparity( screen, DGD3, Param )
end
if Param.StimProtocols.DriftingGratingFlow
    [ DGF, Param ] = SetDriftingGratingsFlow( screen, DGF, Param )
end
if Param.StimProtocols.DriftingGratingFlow2
    [ DGF2, Param ] = SetDriftingGratingsFlow( screen, DGF2, Param )
end
if Param.StimProtocols.DriftingGratingFlow3
    [ DGF3, Param ] = SetDriftingGratingsFlow( screen, DGF3, Param )
end
if Param.StimProtocols.DriftingGratingInDepth
    [ DGID, Param ] = SetDriftingGratingInDepth( screen, DGID, Param )
end
if Param.StimProtocols.DriftingGratingPhi
    [ DGP, Param ] = SetDriftingGratingsPhi( screen(1), DGP, Param )
end
if Param.StimProtocols.StaticGrating
    [ SG, Param ] = SetDriftingGratingsPhi( screen(1), SG, Param )
end
if Param.StimProtocols.RetinotopicPosition
    win = screen(1);
    [ RP , Param ] = SetRetinotopicPosition( win, RP , Param )
    if Param.Dichoptic

    if Param.Dichoptic && Param.stereoMode==0
        win=screen(2);
    else
        win=screen(1);
    end
    [ RP2, Param ] = SetRetinotopicPosition( win, RP2, Param )
    end
end
if Param.StimProtocols.RetinotopicPositionStripes
    win=screen(1);
    [ RPS , Param ] = SetRetinotopicPositionStripes( win, RPS , Param )
    if Param.Dichoptic && Param.stereoMode > 0;
%         if Param.Dichoptic && Param.stereoMode==0
%             win=screen(1);
% %             win=screen(2);
%         else
%             win=screen(1);
%         end
        [ RPS2, Param ] = SetRetinotopicPositionStripes( win, RPS2, Param )
    end
end
if ismember(1, Param.StimProtocols.RFMapping)
    win=screen(1);
    [ RFM, Param ] = SetRFMapping( win, RFM, Param )

    if show_screen_only==1
        for i=1:length(screen)
            mon_pos = get(0,'MonitorPositions');
            figure; set(gcf,'Position',[mon_pos(2,1),mon_pos(2,2),500,500]);
            imshow(RFM.ScreenImg_merge(:,:,:,i),'InitialMagnification','fit','Border','tight'); %#ok<*NODEF>
        end
        Screen('CloseAll'); Priority(0); ShowCursor;
        return
%     else commandwindow;
    end
end
if ismember(2, Param.StimProtocols.RFMapping)
    if Param.Dichoptic && Param.stereoMode==1
        win=screen(1);
    else
        win=screen(2);
    end
    [ RFM2, Param ] = SetRFMapping( win, RFM2, Param )
end
if Param.StimProtocols.PhiMotionGratings
    [ PMG, Param ] = SetPhiMotionGratings( screen(1), PMG, Param )
end
if Param.StimProtocols.DriftingGratingDisparity_phi
    [ DGDphi, Param ] = SetDriftingGratingsDisparity_phi( screen, DGDphi, Param )
end
if Param.StimProtocols.PhiMotionGratingsControl
    [ PMGc1, Param ] = SetPhiMotionGratingsControl( screen(1), PMGc1, Param, 1 )
    [ PMGc2, Param ] = SetPhiMotionGratingsControl( screen(2), PMGc2, Param, 2 )
end

if Param.StimProtocols.PhiGratings
    [ PG, Param ] = SetPhiGratings( screen(1), PG, Param )
end
if Param.StimProtocols.PhiMotionBars
    SetPhiMotionBars
    figure; imshow(ScreenImg_merge(:,:,:,1),'InitialMagnification','fit');
    if Param.Dichoptic
        figure; imshow(ScreenImg_merge(:,:,:,2),'InitialMagnification','fit');
    end
    if show_screen_only==1
        LoadIdentityClut(screen); Screen('CloseAll'); Priority(0); ShowCursor;
        return
    end
end

if Param.StimProtocols.DotMotion
    [ DM , Param ] = SetDotMotion( DM , Param )
end
if Param.StimProtocols.DotMotion2
    [ DM2, Param ] = SetDotMotion( DM2, Param )
end
if Param.StimProtocols.DotMotion3
    [ DM3, Param ] = SetDotMotion( DM3, Param )
end
if Param.StimProtocols.RandomDotStereogram
    [ RDS , Param ] = SetRandomDotStereogram( RDS , Param )
end
if Param.StimProtocols.RandomDotStereogram2
    [ RDS2, Param ] = SetRandomDotStereogram( RDS2, Param )
end
if Param.StimProtocols.RandomBarStereogram
    [ RBS , Param ] = SetRandomBarStereogram( RBS , Param )
end
if Param.StimProtocols.RandomBarStereogram2
    [ RBS2, Param ] = SetRandomBarStereogram( RBS2, Param )
end
if Param.StimProtocols.NaturalImageStereogram
    [ NIS , Param, img ] = SetNaturalImageStereogram( NIS , Param );
end
if Param.StimProtocols.RandomDotStereoBar
    [ RDSB, Param ] = SetRandomDotStereoBar( RDSB, Param )
end
if Param.StimProtocols.RandomDotStereoEdge
    [ RDSE, Param ] = SetRandomDotStereoEdge( RDSE, Param )
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
        if Param.StimProtocols.RetinotopicPosition
            prestimColor = RP.BackgroundLuminance ;
        end
        if Param.StimProtocols.RFMapping
            prestimColor = RFM.BackgroundLuminance ;
        end
        if Param.StimProtocols.RandomDotStereoBar
            prestimColor = RDSB.BackgroundLuminance ;
        end
        if Param.StimProtocols.RandomDotStereoEdge
            prestimColor = RDSE.BackgroundLuminance ;
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
        [DG,  vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, DG );
		    DG.cnt = DG.cnt+1;
    elseif exist('DG2','var') && Param.stimSeq(s)==DG2.stim_id % 
        [DG2, vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, DG2 );
        DG2.cnt = DG2.cnt+1;
    elseif exist('DG3','var') && Param.stimSeq(s)==DG3.stim_id % 
        [DG3, vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, DG3 );
        DG3.cnt = DG3.cnt+1;
    elseif exist('DG4','var') && Param.stimSeq(s)==DG4.stim_id % 
        [DG4, vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, DG4 );
        DG4.cnt = DG4.cnt+1;
    elseif exist('DG5','var') && Param.stimSeq(s)==DG5.stim_id % 
        [DG5, vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, DG5 );
        DG5.cnt = DG5.cnt+1;
    elseif exist('DG6','var') && Param.stimSeq(s)==DG6.stim_id % 
        [DG6, vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, DG6 );
        DG6.cnt = DG6.cnt+1;
    elseif exist('DGD','var') && Param.stimSeq(s)==DGD.stim_id % 
        [DGD, vbl] = PresentDriftingGratingDisparity( nidaq, screen, vbl, Param, DGD );
        DGD.cnt = DGD.cnt+1;
    elseif exist('DGD2','var') && Param.stimSeq(s)==DGD2.stim_id % 
        [DGD2, vbl] = PresentDriftingGratingDisparity( nidaq, screen, vbl, Param, DGD2 );
        DGD2.cnt = DGD2.cnt+1;
    elseif exist('DGD3','var') && Param.stimSeq(s)==DGD3.stim_id % 
        [DGD3, vbl] = PresentDriftingGratingDisparity( nidaq, screen, vbl, Param, DGD3 );
        DGD3.cnt = DGD3.cnt+1;
    elseif exist('DGF','var') && Param.stimSeq(s)==DGF.stim_id % 
        [DGF, vbl] = PresentDriftingGratingFlow( nidaq, screen, vbl, Param, DGF );
        DGF.cnt = DGF.cnt+1;
    elseif exist('DGF2','var') && Param.stimSeq(s)==DGF2.stim_id % 
        [DGF2, vbl] = PresentDriftingGratingFlow( nidaq, screen, vbl, Param, DGF2 );
        DGF2.cnt = DGF2.cnt+1;
    elseif exist('DGF3','var') && Param.stimSeq(s)==DGF3.stim_id % 
        [DGF3, vbl] = PresentDriftingGratingFlow( nidaq, screen, vbl, Param, DGF3 );
        DGF3.cnt = DGF3.cnt+1;
    elseif exist('DGID','var') && Param.stimSeq(s)==DGID.stim_id % 
        [DGID, vbl] = PresentDriftingGratingInDepth( nidaq, screen, vbl, Param, DGID );
        DGID.cnt = DGID.cnt+1;
    elseif exist('DGP','var') && Param.stimSeq(s)==DGP.stim_id % 
        [DGP, vbl] = PresentDriftingGratingPhi( nidaq, screen, vbl, Param, DGP );
        DGP.cnt = DGP.cnt+1;
    elseif exist('SG','var') && Param.stimSeq(s)==SG.stim_id % 
        [SG, vbl] = PresentDriftingGratingPhi( nidaq, screen, vbl, Param, SG );
        SG.cnt = SG.cnt+1;
	elseif exist('RP','var') && Param.stimSeq(s)==RP.stim_id %
        if Param.Dichoptic && Param.stereoMode == 0
            Screen('FillRect', screen(2), 127);
            Screen('Flip', screen(2));
        end
		[RP, vbl] = PresentRetinotopicPosition( nidaq, screen(1), vbl, Param, RP );
		RP.cnt = RP.cnt+1;
	elseif exist('RP2','var') && Param.stimSeq(s)==RP2.stim_id %
        if Param.Dichoptic && Param.stereoMode == 0
            Screen('FillRect', screen(1), 127);
            Screen('Flip', screen(1));
            [RP2, vbl] = PresentRetinotopicPosition( nidaq, screen(2), vbl, Param, RP2 );
        else
            [RP2, vbl] = PresentRetinotopicPosition( nidaq, screen(1), vbl, Param, RP2 );
        end
		RP2.cnt = RP2.cnt+1;
	elseif exist('RPS','var') && Param.stimSeq(s)==RPS.stim_id %
        if Param.Dichoptic && Param.stereoMode == 0
%             Screen('FillRect', screen(2), 127);
%             Screen('Flip', screen(2));
        end
		[RPS, vbl] = PresentRetinotopicPositionStripes( nidaq, win, vbl, Param, RPS );
		RPS.cnt = RPS.cnt+1;
	elseif exist('RPS2','var') && Param.stimSeq(s)==RPS2.stim_id %
        if Param.Dichoptic && Param.stereoMode == 0
%             Screen('FillRect', screen(1), 127);
%             Screen('Flip', screen(1));
        end
		[RPS2, vbl] = PresentRetinotopicPositionStripes( nidaq, win, vbl, Param, RPS2 );
		RPS2.cnt = RPS2.cnt+1;
    elseif exist('RFM','var') && Param.stimSeq(s)==RFM.stim_id %
        [RFM, vbl] = PresentRFMapping( nidaq, screen, vbl, Param, RFM );
        RFM.cnt = RFM.cnt+1;
    elseif exist('RFM2','var') && Param.stimSeq(s)==RFM2.stim_id %
        [RFM2, vbl] = PresentRFMapping( nidaq, screen, vbl, Param, RFM2 );
        RFM2.cnt = RFM2.cnt+1;
    elseif exist('PMG','var') && Param.stimSeq(s)==PMG.stim_id %
		[PMG, vbl] = PresentPhiMotionGrating( nidaq, screen, vbl, Param, PMG );
		PMG.cnt = PMG.cnt+1;
    elseif exist('DGDphi','var') && Param.stimSeq(s)==DGDphi.stim_id % 
        [DGDphi, vbl] = PresentDriftingGratingDisparity_phi( nidaq, screen, vbl, Param, DGDphi );
        DGDphi.cnt = DGDphi.cnt+1;
    elseif exist('PMGc1','var') && Param.stimSeq(s)==PMGc1.stim_id %
		[PMGc1, vbl] = PresentPhiMotionGratingControl( nidaq, screen, vbl, Param, PMGc1, 1 );
		PMGc1.cnt = PMGc1.cnt+1;
    elseif exist('PMGc2','var') && Param.stimSeq(s)==PMGc2.stim_id %
		[PMGc2, vbl] = PresentPhiMotionGratingControl( nidaq, screen, vbl, Param, PMGc2, 2 );
		PMGc2.cnt = PMGc2.cnt+1;
	elseif exist('PG','var') && Param.stimSeq(s)==PG.stim_id %
		[PG, vbl] = PresentPhiGratings( nidaq, screen, vbl, Param, PG );
		PG.cnt = PG.cnt+1;
    elseif exist('PMB','var') && Param.stimSeq(s)==PMB.stim_id %
		[PMB, vbl] = PresentPhiMotionBars( nidaq, screen, vbl, Param, PMB );
		PMB.cnt = PMB.cnt+1;
    elseif exist('DM','var') && Param.stimSeq(s)==DM.stim_id % 
        [DM,  vbl] = PresentDotMotion( nidaq, screen, vbl, Param, DM );
		DM.cnt = DM.cnt+1;
    elseif exist('DM2','var') && Param.stimSeq(s)==DM2.stim_id % 
        [DM2,  vbl] = PresentDotMotion( nidaq, screen, vbl, Param, DM2 );
		DM2.cnt = DM2.cnt+1;
    elseif exist('DM3','var') && Param.stimSeq(s)==DM3.stim_id % 
        [DM3,  vbl] = PresentDotMotion( nidaq, screen, vbl, Param, DM3 );
		DM3.cnt = DM3.cnt+1;
    elseif exist('RDS','var') && Param.stimSeq(s)==RDS.stim_id % 
        [RDS,  vbl] = PresentRandomDotStereogram( nidaq, screen, vbl, Param, RDS );
		RDS.cnt = RDS.cnt+1;
    elseif exist('RDS2','var') && Param.stimSeq(s)==RDS2.stim_id % 
        [RDS2,  vbl] = PresentRandomDotStereogram( nidaq, screen, vbl, Param, RDS2 );
		RDS2.cnt = RDS2.cnt+1;
    elseif exist('RBS','var') && Param.stimSeq(s)==RBS.stim_id % 
        [RBS,  vbl] = PresentRandomBarStereogram( nidaq, screen, vbl, Param, RBS );
		RBS.cnt = RBS.cnt+1;
    elseif exist('RBS2','var') && Param.stimSeq(s)==RBS2.stim_id % 
        [RBS2,  vbl] = PresentRandomBarStereogram( nidaq, screen, vbl, Param, RBS2 );
		RBS2.cnt = RBS2.cnt+1;
    elseif exist('NIS','var') && Param.stimSeq(s)==NIS.stim_id % 
        [NIS,  vbl] = PresentNaturalImageStereogram( nidaq, screen, vbl, Param, NIS, img );
		NIS.cnt = NIS.cnt+1;
    elseif exist('RDSB','var') && Param.stimSeq(s)==RDSB.stim_id % 
        [RDSB,  vbl] = PresentRandomDotStereoBar( nidaq, screen, vbl, Param, RDSB );
		RDSB.cnt = RDSB.cnt+1;
    elseif exist('RDSE','var') && Param.stimSeq(s)==RDSE.stim_id % 
        [RDSE,  vbl] = PresentRandomDotStereoEdge( nidaq, screen, vbl, Param, RDSE );
		RDSE.cnt = RDSE.cnt+1;
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

SaveParameters_script;
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
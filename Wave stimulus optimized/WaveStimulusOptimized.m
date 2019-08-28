function WaveStimulusOptimized

StimSettings.ScreenDist = 18;
Settings.NumTrials = 20;
StimSettings.WaveType = 'FlickerShifty';
StimSettings.BarSize = 20;
StimSettings.BarSpeed = 18; 
StimSettings.prestimInterval = 20;
StimSettings.LoopDelay = 3; % sec
StimSettings.GratShape = 'Square';
StimSettings.SpatFreq = 0.04;
StimSettings.TempFreq = 2;
StimSettings.FlickerFreq = 6;
StimSettings.MinInt = 0;
StimSettings.MaxInt = 255;
StimSettings.BGint = 127;
StimSettings.ScrResX = 1600;
StimSettings.ScrResY =  900;
StimSettings.ScrSizeX = 44;
StimSettings.ScrSizeY = 25;
StimSettings.PresAreaX = StimSettings.ScrResX;
StimSettings.PresAreaY = StimSettings.ScrResY;
Settings.DipslayBlueGreen = 0;
StimSettings.XYsimult = 0;
Settings.useArduino = 1;
loadMyGammaTable = 1;
CorrectCurvature = 1;
Settings.test = 0;

%% Load Gamma table
if loadMyGammaTable
    try
        % get calibrated gammaTable:
        if     strfind(eval('computer'),'WIN')
            load('C:\Users\lachioma\Dropbox\Code\alessandro\screen_calibration\2016-05-04\MyGammaTable.mat');
        elseif strfind(eval('computer'),'LNX')
            load('/home/alessandro/Dropbox/Code/alessandro/screen_calibration/2016-05-04/MyGammaTable.mat');
        end
        Settings.GammaCorrectionTable = gammaTable;
    catch
        disp('WARNING: No gamma correction implemented!!');
        Settings.GammaCorrectionTable = [];
    end
end
%% Set save file path
if     strfind(eval('computer'),'WIN')
    save_dir=['C:\Code\alessandro\stim_data\' datestr(now, 'yyyy_mm_dd')];
elseif strfind(eval('computer'),'LNX')
    save_dir=['/home/alessandro/Dropbox/Code/alessandro/stim_data/' datestr(now, 'yyyy_mm_dd')];
end
if ~isdir(save_dir), mkdir(save_dir), end;
datetime_suffix = datestr(now, 'yyyy_mm_dd_HH_MM_SS') ;
stim_fname = [save_dir filesep 'InfoStim-' datetime_suffix '-WaveStimulusOptimized'];
Settings.stim_fname = stim_fname;
%% Set Arduino
if Settings.useArduino == 1 && Settings.test == 0
%     nidaq=[];
    Settings.ArdPinNr = 9;
    global Ard %#ok<TLEV>
    if isempty(Ard)
        Ard = arduino('/dev/ttyUSB0');
        Ard.pinMode(Settings.ArdPinNr,'output');
    else
        disp('Arduino already connected !');
    end
    Ard.digitalWrite(Settings.ArdPinNr,0); % set to 0V
else
    Settings.ArdPinNr = [];
end
%% Set curvature correctionsca
if CorrectCurvature
    warptype = 'CurvatureCorrection';
    rotationAngle = 0;
    inSize = [StimSettings.ScrResX,StimSettings.ScrResY];
    inOffset = [0,0];
    outSize = [StimSettings.ScrResX,StimSettings.ScrResY];
    outOffset = [0,0];
    ScrSz = [StimSettings.ScrSizeX StimSettings.ScrSizeY];
    ScrDist = StimSettings.ScreenDist;
    Wflat = 20;
    R = 12;
    SaveName = [save_dir filesep 'RM_CurvCorrCalib-' datetime_suffix '.mat'];
    save(SaveName,'warptype','rotationAngle','inSize','inOffset',...
        'outSize','outOffset','ScrSz', 'ScrDist', 'Wflat', 'R', '-mat', '-V6');
    Settings.CurvCorrectFile = SaveName;
else
    Settings.CurvCorrectFile = [];
end

%% Run Wave
RunWave(Settings,StimSettings);



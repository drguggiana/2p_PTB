
function [stimtype, vbl] = PresentLooming( nidaq, screen, vbl, Param, stimtype, tar_eye )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
tic;

if     Param.Dichoptic && Param.stereoMode==0 
%     if     strcmp(inputname(5),'DG') || strcmp(inputname(5),'DG3')
        win=screen(tar_eye);
%     elseif strcmp(inputname(5),'DG2')|| strcmp(inputname(5),'DG4')
%         win=screen(2);
%     end
    winStereo = 0; % for Photodiode
    Screen('FillRect', screen(1), stimtype.BackgroundLuminance);
    Screen('FillRect', screen(2), stimtype.BackgroundLuminance);
    vbl = Screen('Flip', screen(1));
    vbl = Screen('Flip', screen(2));
% elseif Param.Dichoptic && ismember(Param.stereoMode,[1 4])
%     if     strcmp(inputname(5),'DG')  || ...
%            strcmp(inputname(5),'DG3') || ...
%            strcmp(inputname(5),'DG5')
%         winStereo = 0;
%     elseif strcmp(inputname(5),'DG2') || ...
%            strcmp(inputname(5),'DG4') || ...
%            strcmp(inputname(5),'DG6')
%         winStereo = 1;
%     end
%     win=screen(1);
%     Screen('SelectStereoDrawBuffer', win, 0);
%     Screen('FillRect', win, stimtype.BackgroundLuminance);
%     Screen('SelectStereoDrawBuffer', win, 1);
%     Screen('FillRect', win, stimtype.BackgroundLuminance);
%     vbl = Screen('Flip', win);
else
    win=screen(1);
    winStereo = 0;
    vbl = Screen('Flip', win);
end

% if stimtype.drawMask == 0
% 	Screen('BlendFunction', win, GL_ONE, GL_ZERO);
% else % if gaussian transparency mask:
% 	% Enable alpha blending for combination of the gaussian aperture
% 	% with the drifting grating:
% 	Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% end
% Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
% % Disable alpha-blending, restrict following drawing to alpha channel:
% Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]);

%define the rectangle that will bind the circle, using its max size
baseRect = [0 0 1 1];

[xCenter, yCenter] = RectCenter(stimtype.gratingRect);
% Center the rectangle on the centre of the screen
centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
%define the max Diameter
%max_diam = max(screenYpixels);
% Sync us and get a time stamp
vbl = Screen('Flip', win);

%define the coordinates of the circle center
circle_center_x = 0.5*stimtype.gratingsize(1);
circle_center_y = 0.5*stimtype.gratingsize(2);

% %define the final size of the circle
% max_size = 0.33*screenXpixels;

% %define the scale factor in pixels/frame
% scaleFactor = max_size/n_frames;


test=Param.test;
tempo=[0 0];
% seqangles = stimtype.seqangles';
% angle = seqangles(stimtype.cnt)
% if ~isfield(stimtype,'phase') || isempty(stimtype.phase)
%     phase = round((359-0)*rand(1) + 0);
% else
%     phase = stimtype.phase;
% end

%transpose the times and speeds arrays
seqtimes = stimtype.seqtimes';
seqspeeds = stimtype.seqspeeds';
seqspeedorder = stimtype.paramorder(:,:,1)';

vblendtime = vbl + seqtimes(stimtype.cnt) - Param.ifi;

%define the scale factor for the expanding circle
scaleFactor = seqspeeds(stimtype.cnt);
%get the current color
stim_color = stimtype.seqcolors(stimtype.cnt,:);
%also print the current expansion speed
fprintf(strcat('Current exp speed:',...
    num2str(stimtype.expansion_speeds(seqspeedorder(stimtype.cnt))),...
    '/Current color:',num2str(stimtype.seqcolors(stimtype.cnt,:)),...
    '\r\n'))
%         for i = 1 : frame_stimulus_dur
frameNr = 0;

while vbl < vblendtime
    frameNr = frameNr + 1;
	% Update some grating animation parameters:
% 	% Increment phase by 1 degree:
% 	phase = phase + stimtype.phaseincrement;
	% Draw the grating, centered on the screen, with given rotation 'angle',
	% sine grating 'phase' shift and amplitude, rotating via set
	% 'rotateMode'. Note that we pad the last argument with a 4th
	% component, which is 0. This is required, as this argument must be a
	% vector with a number of components that is an integral multiple of 4,
	% i.e. in our case it must have 4 components:
%     if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
%         Screen('SelectStereoDrawBuffer', win, winStereo);
%     end
% 	Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingRect, angle, [], [], [], [], Param.rotateMode, [phase, stimtype.freq_cppx, stimtype.amplitude, 0]);
    
%variable intensity circle
% Screen('FillOval', window, [1 1 1]*(n_frames-frames)/n_frames, centeredRect); 
    Screen('FillOval', win,stim_color, centeredRect);
     
    %update the size of the circle
    baseRect(3:4) = baseRect(3:4) + scaleFactor;
    
    % Center the rectangle on the centre of the screen
    centeredRect = CenterRectOnPointd(baseRect, circle_center_x, circle_center_y);
    
% 	if stimtype.drawMask == 1
%         Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% 		% Draw mask over grating:
% 		Screen('DrawTexture', win, stimtype.masktex{1}, [], stimtype.gratingRect);
%         % Disable alpha-blending:
%         Screen('Blendfunction', win, GL_ONE, GL_ZERO);
% 	end;
    
    if Param.usePhotodiode
        if ismember(Param.stereoMode,[1 4])
            Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
        end
        Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
    end
    
    Screen('DrawingFinished', win);

    if frameNr == 1
        StimSignal( Param, test, nidaq, stimtype.stim_id )
    end
    
    % Show it at next retrace:
	vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);
    
    Exit_If_Shift(test,nidaq);
end
tempo(1) = toc;

StimSignal( Param, test, nidaq, 0 )

% Post-stim:
tic
Screen('FillRect', win, stimtype.BackgroundLuminance);
while  tempo(2) < stimtype.poststim_time
	tempo(2) = toc;
	Exit_If_Shift(test,nidaq);
    vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);  
end
stimtype.tempo{stimtype.cnt}=tempo;
	
if exist('scriptName','var') && stimtype.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end

function [stimtype, vbl] = PresentDriftingGrating( nidaq, screen, vbl, Param, stimtype, tar_eye )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
tic;

if     Param.Dichoptic && Param.stereoMode==0 
%     if     strcmp(inputname(5),'DG') || strcmp(inputname(5),'DG3')
%         win=screen(1);
%     elseif strcmp(inputname(5),'DG2')|| strcmp(inputname(5),'DG4')
%         win=screen(2);
%     end
    win = screen(tar_eye);
    winStereo = 0; % for Photodiode
    Screen('FillRect', screen(1), stimtype.BackgroundLuminance);
    Screen('FillRect', screen(2), stimtype.BackgroundLuminance);
    vbl = Screen('Flip', screen(1));
    vbl = Screen('Flip', screen(2));
elseif Param.Dichoptic && ismember(Param.stereoMode,[1 4])
    if     strcmp(inputname(5),'DG')  || ...
           strcmp(inputname(5),'DG3') || ...
           strcmp(inputname(5),'DG5')
        winStereo = 0;
    elseif strcmp(inputname(5),'DG2') || ...
           strcmp(inputname(5),'DG4') || ...
           strcmp(inputname(5),'DG6')
        winStereo = 1;
    end
    win=screen(1);
    Screen('SelectStereoDrawBuffer', win, 0);
    Screen('FillRect', win, stimtype.BackgroundLuminance);
    Screen('SelectStereoDrawBuffer', win, 1);
    Screen('FillRect', win, stimtype.BackgroundLuminance);
    vbl = Screen('Flip', win);
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


test=Param.test;
tempo=[0 0];
seqangles = stimtype.seqangles;
angle = seqangles(stimtype.cnt);
if ~isfield(stimtype,'phase') || isempty(stimtype.phase)
    phase = round((359-0)*rand(1) + 0);
else
    phase = stimtype.phase;
end
vblendtime = vbl + stimtype.stimulus_time - Param.ifi;
%         for i = 1 : frame_stimulus_dur
frameNr = 0;

%report the current stimulus
fprintf(strcat('Angle:',num2str(angle),...
    '/SF:',num2str(stimtype.seq_spatial(stimtype.cnt)),...
    '/TF:',num2str(stimtype.seq_temporal(stimtype.cnt)),'\r\n'))

while vbl < vblendtime
    frameNr = frameNr + 1;
	% Update some grating animation parameters:
	% Increment phase by 1 degree:
	phase = phase + stimtype.seq_temporal(stimtype.cnt);
	% Draw the grating, centered on the screen, with given rotation 'angle',
	% sine grating 'phase' shift and amplitude, rotating via set
	% 'rotateMode'. Note that we pad the last argument with a 4th
	% component, which is 0. This is required, as this argument must be a
	% vector with a number of components that is an integral multiple of 4,
	% i.e. in our case it must have 4 components:
    if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
        Screen('SelectStereoDrawBuffer', win, winStereo);
    end
	Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingRect, angle, [], [], [], [], Param.rotateMode,...
        [phase, stimtype.seq_spatial(stimtype.cnt), stimtype.amplitude, 0]);
    
	if stimtype.drawMask == 1
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		% Draw mask over grating:
		Screen('DrawTexture', win, stimtype.masktex{1}, [], stimtype.gratingRect);
        % Disable alpha-blending:
        Screen('Blendfunction', win, GL_ONE, GL_ZERO);
    end
    
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
while  tempo(2) < stimtype.poststim_time;
	tempo(2) = toc;
	Exit_If_Shift(test,nidaq);
    vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);  
end
stimtype.tempo{stimtype.cnt}=tempo;
	
    if exist('scriptName','var') && stimtype.cnt==1
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
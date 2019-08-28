
function [stimtype, vbl] = PresentDriftingGratingFlow( nidaq, screen, vbl, Param, stimtype )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
tic;

win=screen(1);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);

Screen('SelectStereoDrawBuffer', win, 0);
Screen('FillRect', win, stimtype.BackgroundLuminance);
Screen('SelectStereoDrawBuffer', win, 1);
Screen('FillRect', win, stimtype.BackgroundLuminance);
vbl = Screen('Flip', win);


test=Param.test;
tempo=[0 0];
cnt=stimtype.cnt;
stimseq=stimtype.stimseq';
seqangles = stimtype.seqangles';
seqphases = stimtype.seqphases';
angle0 = seqangles(stimseq(cnt));
if angle0 >= 180
    angle1 = angle0 - 180 + stimtype.offset_rot_deg(1);
    angle2 = angle0 - 180 + stimtype.offset_rot_deg(2);
    phaseincrement = - stimtype.phaseincrement;
else
    angle1 = angle0 + stimtype.offset_rot_deg(1);
    angle2 = angle0 + stimtype.offset_rot_deg(2);
    phaseincrement = stimtype.phaseincrement;
end

phase1 = stimtype.InitialPhase(cnt);
phase2 = phase1;
deltaphase = stimtype.IOPhaseDifferences(seqphases(stimseq(cnt)));
disp(['Angle: ' num2str(seqangles(stimseq(cnt))) ', Phase1: ' num2str(phase1) ', DeltaPhase: ' num2str(deltaphase)]);

vblendtime = vbl + stimtype.stimulus_time - Param.ifi;
%         for i = 1 : frame_stimulus_dur
frameNr = 0;

while vbl < vblendtime
    frameNr = frameNr + 1;
	% Update some grating animation parameters:
	% Increment phase by 1 degree:
% 	stimtype.phase = stimtype.phase + stimtype.phaseincrement;
    phase1 = phase1 + phaseincrement;
    phase2 = phase2 - phaseincrement;
    
	% Draw the grating, centered on the screen, with given rotation 'angle',
	% sine grating 'phase' shift and amplitude, rotating via set
	% 'rotateMode'. Note that we pad the last argument with a 4th
	% component, which is 0. This is required, as this argument must be a
	% vector with a number of components that is an integral multiple of 4,
	% i.e. in our case it must have 4 components:
    
    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);
    eye = 1;
	Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingRect(eye,:), angle1, [], [], [], [], Param.rotateMode, [phase1           , stimtype.freq_cppx(1), stimtype.amplitude, 0]);
    if stimtype.drawMask == 1
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		% Draw mask over grating:
		Screen('DrawTexture', win, stimtype.masktex{eye}, [], stimtype.gratingRect(eye,:));
        % Disable alpha-blending:
        Screen('Blendfunction', win, GL_ONE, GL_ZERO);
	end;
	if Param.usePhotodiode
        Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
		Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
	end
    Screen('SelectStereoDrawBuffer', win, 1);
    eye = 1;
	Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingRect(eye,:), angle2, [], [], [], [], Param.rotateMode, [phase2+deltaphase, stimtype.freq_cppx(2), stimtype.amplitude, 0]);
    if stimtype.drawMask == 1
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		% Draw mask over grating:
		Screen('DrawTexture', win, stimtype.masktex{eye}, [], stimtype.gratingRect(eye,:));
        % Disable alpha-blending:
        Screen('Blendfunction', win, GL_ONE, GL_ZERO);
    end;
    Screen('DrawingFinished', win);
    
	if frameNr == 1 %#ok<ALIGN>
		StimSignal( Param, test, nidaq, stimtype.stim_id )
    end
    
    % Show it at next retrace:
	vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);
%     vbl = Screen('Flip', win );
    Exit_If_Shift(test,nidaq);
end
tempo(1) = toc;
StimSignal( Param, test, nidaq, 0 )
tic
Screen('FillRect', win, stimtype.BackgroundLuminance);

while  tempo(2) < stimtype.poststim_time;
	tempo(2) = toc;
	Exit_If_Shift(test,nidaq);
    vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);  
end
stimtype.tempo{cnt}=tempo;
	
    if exist('scriptName','var') && stimtype.cnt==1
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
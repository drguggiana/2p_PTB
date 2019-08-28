
function [stimtype, vbl] = PresentDriftingGratingInDepth( nidaq, screen, vbl, Param, stimtype )
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
stimseq       = stimtype.stimseq';
seqangles     = stimtype.seqangles';
seqconditions = stimtype.seqconditions';
angle0 = seqangles(stimseq(cnt));

angle1 = angle0 + stimtype.offset_rot_deg(1);
angle2 = angle0 + stimtype.offset_rot_deg(2);
phaseincrement1 = stimtype.phaseincrement(seqconditions(stimseq(cnt)),1);
phaseincrement2 = stimtype.phaseincrement(seqconditions(stimseq(cnt)),2);

InitialPhase1 = stimtype.InitialPhase(:,:,1);
InitialPhase2 = stimtype.InitialPhase(:,:,2);
phase1 = InitialPhase1(cnt);
phase2 = InitialPhase2(cnt);

% % deltaphase = stimtype.IOPhaseDifferences(seqphases(stimseq(cnt)));
% % disp(['Angle: ' num2str(seqangles(stimseq(cnt))) ', Phase1: ' num2str(phase1) ', DeltaPhase: ' num2str(deltaphase)]);
disp(['Angle: ' num2str(angle0) ', speed1: ' num2str(stimtype.speedCombinations(seqconditions(stimseq(cnt)),1)) ', speed2: ' num2str(stimtype.speedCombinations(seqconditions(stimseq(cnt)),2))]);

% %initial phases corrected for RF centering:
% phase1 = phase1 + abs(cos(rad(angle1)))*stimtype.offset_phase1_deg(1) + abs(sin(rad(angle1)))*stimtype.offset_phase1_deg(2);
% phase2 = phase2 + abs(cos(rad(angle2)))*stimtype.offset_phase2_deg(1) + abs(sin(rad(angle2)))*stimtype.offset_phase2_deg(2);

phase1running = phase1 - phaseincrement1;
phase2running = phase2 - phaseincrement2;

vblendtime = vbl + stimtype.stimulus_time - Param.ifi;
%         for i = 1 : frame_stimulus_dur
frameNr = 0;

while vbl < vblendtime
    frameNr = frameNr + 1;
	% Update grating phase:
%     phase1 = phase1 + phaseincrement;
    phase1running = phase1running + phaseincrement1;
    phase2running = phase2running + phaseincrement2;
    
	% Draw the grating, centered on the screen, with given rotation 'angle',
	% sine grating 'phase' shift and amplitude, rotating via set
	% 'rotateMode'. Note that we pad the last argument with a 4th
	% component, which is 0. This is required, as this argument must be a
	% vector with a number of components that is an integral multiple of 4,
	% i.e. in our case it must have 4 components:
    
    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);
    eye = 1;
	Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingRect(eye,:), angle1, [], [], [], [], Param.rotateMode, [phase1running, stimtype.freq_cppx(1), stimtype.amplitude, 0]);
    if stimtype.drawMask == 1
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		% Draw mask over grating:
		Screen('DrawTexture', win, stimtype.masktex{eye}, [], stimtype.gratingRect(eye,:));
%         Screen('DrawTexture', win, stimtype.masktex{eye}, [], []);
        % Disable alpha-blending:
        Screen('Blendfunction', win, GL_ONE, GL_ZERO);
	end;
	if Param.usePhotodiode
        Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
		Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
	end
    Screen('SelectStereoDrawBuffer', win, 1);
    eye = 1;
	Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingRect(eye,:), angle2, [], [], [], [], Param.rotateMode, [phase2running, stimtype.freq_cppx(2), stimtype.amplitude, 0]);
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
    
        function alpharad = rad(alphadeg)
            alpharad=alphadeg*pi/180;
        end

end

function [stimtype, vbl] = PresentDriftingGratingDisparity_phi( nidaq, screen, vbl, Param, stimtype )
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
tempo=[];
cnt=stimtype.cnt;
stimseq=stimtype.stimseq';
seqangles = stimtype.seqangles';
seqphases = stimtype.seqphases';
angle0 = seqangles(stimseq(cnt));
anglecart = mod(angle0-90, 360);
if angle0 >= 180
    angle1 = angle0 - 180 + stimtype.offset_rot_deg(1);
    angle2 = angle0 - 180 + stimtype.offset_rot_deg(2);
    phaseincrement = - stimtype.phaseincrement;
else
    angle1 = angle0 + stimtype.offset_rot_deg(1);
    angle2 = angle0 + stimtype.offset_rot_deg(2);
    phaseincrement = stimtype.phaseincrement;
end

InitialPhase = stimtype.InitialPhase(stimseq(cnt));
% this is to take into account the following add of phaseincrement, so that
% the first phase to be displayed will really be InitialPhase(stimseq(cnt)):
phase1 = InitialPhase - phaseincrement;
phase2 = phase1;
deltaphase = stimtype.IOPhaseDifferences(seqphases(stimseq(cnt)));
fprintf('  * Angle=%3.0f, (cartesian %3.0f); Phase1=%3.0f, DeltaPhase=%3.0f \n',angle0,anglecart,InitialPhase,deltaphase);

phase1 = phase1 + abs(cos(rad(angle1)))*stimtype.offset_phase1_deg(1) + abs(sin(rad(angle1)))*stimtype.offset_phase1_deg(2);
phase2 = phase2 + abs(cos(rad(angle2)))*stimtype.offset_phase2_deg(1) + abs(sin(rad(angle2)))*stimtype.offset_phase2_deg(2);

vblendtime = vbl + stimtype.stimulus_time - Param.ifi;
go_on = 1;

frameNr = 0;
patchPresentationNr = 0;

while go_on == 1 
    vblendpatch = vbl + stimtype.patch_time - Param.ifi;
    patchPresentationNr = patchPresentationNr + 1;
    % Update grating phase:
    phase1 = phase1 + phaseincrement;
    phase2 = phase2 + phaseincrement;
    
    while vbl < vblendpatch
        frameNr = frameNr + 1;

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

        if frameNr == 1
            StimSignal( Param, test, nidaq, stimtype.stim_id )
        end

        % Show it at next retrace:
        vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);
    %     vbl = Screen('Flip', win );
        if vbl > vblendtime
            go_on=0;
            break
        end
        Exit_If_Shift(test,nidaq);
    end
    tempo = [tempo toc]; %#ok<*AGROW>
    
end
StimSignal( Param, test, nidaq, 0 )
tic
Screen('FillRect', win, stimtype.BackgroundLuminance);

while  vbl < vblendtime + stimtype.poststim_time
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
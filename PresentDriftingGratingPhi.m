
function [stimtype, vbl] = PresentDriftingGratingPhi( nidaq, screen, vbl, Param, stimtype )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
tic;
if     Param.Dichoptic && strcmp(inputname(5),'DG')
    win=screen(1);
    if Param.StimProtocols.RetinotopicPosition || Param.StimProtocols.PhiMotion
        Screen('FillRect', screen(2), 127);
    else
        Screen('FillRect', screen(2), 0);
    end
    Screen('Flip', screen(2));
elseif Param.Dichoptic && strcmp(inputname(5),'DG2')
    win=screen(2);
    if Param.StimProtocols.RetinotopicPosition || Param.StimProtocols.PhiMotion
        Screen('FillRect', screen(1), 127);
    else
        Screen('FillRect', screen(1), 0);
    end
    Screen('Flip', screen(1));
else
    win=screen(1);
end
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
test=Param.test;
tempo=[0 0 0];
cnt=stimtype.cnt;
stimseq=stimtype.stimseq';
seqangles = stimtype.seqangles';
seqphases = stimtype.seqphases';
angle = seqangles(stimseq(cnt));
phase = stimtype.phase(seqphases(stimseq(cnt)));
if angle >= 180
    angle = angle - 180;
    phaseincrement = - stimtype.phaseincrement;
else
    phaseincrement = + stimtype.phaseincrement;
end
% vblendtime = vbl + 2*stimtype.patch_time - Param.ifi;  %Ale: I put -ifi to obtain an exact stimulus time, I dont know why exactly.
%         for i = 1 : frame_stimulus_dur
for p = 1:2
    if p == 2
        phase = phase + phaseincrement;
    end
    for f = 1 : stimtype.patch_time_fr

        % Draw the grating, centered on the screen, with given rotation 'angle',
        % sine grating 'phase' shift and amplitude, rotating via set
        % 'rotateMode'. Note that we pad the last argument with a 4th
        % component, which is 0. This is required, as this argument must be a
        % vector with a number of components that is an integral multiple of 4,
        % i.e. in our case it must have 4 components:
        Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase, stimtype.freq_cppx, stimtype.amplitude, 0]);
		if Param.usePhotodiode
			Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
        end
        
        if p==1 && f==1
            StimSignal( Param, test, nidaq, stimtype.stim_id )
        end
        
        % Show it at next retrace:
        vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);
        Exit_If_Shift(test,nidaq);
    end
    tempo(p) = toc;
end
StimSignal( Param, test, nidaq, 0 )
tic
Screen('FillRect', win, stimtype.BackgroundLuminance);

while  tempo(3) < stimtype.poststim_time;
	tempo(3) = toc;
	Exit_If_Shift(test,nidaq);
    vbl = Screen('Flip', win);  
end
stimtype.tempo{cnt}=tempo;


if exist('scriptName','var') && stimtype.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end

function [stimtype, vbl] = PresentPhiGratings( nidaq, screen, vbl, Param, stimtype )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

for i=1:length(screen)
    Screen('BlendFunction', screen(i), GL_ONE, GL_ZERO);
    Screen('FillRect', screen(i), stimtype.BackgroundLuminance);
    vbl = Screen('Flip', screen(i)); 
end
cnt=stimtype.cnt;
stimseq=stimtype.stimseq';
test=Param.test;
stimtype.tempo(cnt).toc_poststim = 0;
% angle = stimtype.seqangles(stimseq(cnt));
tempo=[];
frameNr=0;

    for s = 1:2
        f = stimtype.screensequence(stimseq(cnt),s);
%         for f = stimtype.screensequence(stimseq(cnt),:) % [1,2] or [2,1]
        if f==0; continue; end;
        tic
        
        angle = stimtype.seqangles(stimseq(cnt)) + stimtype.offset_rot_deg(f);
        phase = stimtype.phases(stimtype.seqphases(stimseq(cnt),s));

        vblendpatch = vbl + stimtype.patch_time - Param.ifi;
        while vbl < vblendpatch
            frameNr = frameNr + 1;
            Screen('DrawTexture', screen(f), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase, stimtype.freq_cppx(f), stimtype.amplitude, 0]);
%                     Screen('FillRect', screen(f), stimtype.BackgroundLuminance );
%                     vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
%                     vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi );
			if Param.usePhotodiode
				Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
            end
            if frameNr==1
                StimSignal( Param, test, nidaq, stimtype.stim_id )
            end
            vbl  = Screen('Flip', screen(f) );
%                     vbl2 = Screen('Flip', screen(2) );
%                     vbl=vbl2;

            Exit_If_Shift(test,nidaq);  
        end
        tempo = [tempo toc]; %#ok<*AGROW>
        
        tic_interpatch = tic;
        
        % If the stimulus is not the same phase presented in the same
        % screen, fill the screen that has just been presented and flip it
        % (otherwise the texture remains on the screen if it's not
        % the same screen to be presented in s=2).
        % In addition, present interpatch_time (only when >0.017 because
        % the previous filling already takes one refresh time)
        if diff( stimtype.screensequence(stimseq(cnt),:) )==0 &&...
                diff( stimtype.seqphases(stimseq(cnt),:) )==0
        else
            Screen('FillRect', screen(f), stimtype.BackgroundLuminance);
            vbl = Screen('Flip', screen(f)); 
            % interpatch
            if s==1 && stimtype.interpatch_time > 0.017 %&& diff( stimtype.screensequence(stimseq(cnt),:) )~=0
                Screen('FillRect', screen(1), stimtype.BackgroundLuminance);
				if Param.usePhotodiode
					Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
				end
                Screen('FillRect', screen(2), stimtype.BackgroundLuminance);
                for fr = 1 : stimtype.frame_interpatch
                    Screen('Flip', screen(1) );
                    vbl = Screen('Flip', screen(2) );
                end
            end
        end
        tempo = [tempo toc(tic_interpatch)];
        
    end


    %%% post-stim blank screen
    tic_poststim = tic;
    
    StimSignal( Param, test, nidaq, 0 )

    Screen('FillRect', screen(1),stimtype.BackgroundLuminance);
    Screen('FillRect', screen(2),stimtype.BackgroundLuminance);
    while stimtype.tempo(cnt).toc_poststim < stimtype.poststim_time
        stimtype.tempo(cnt).toc_poststim = toc(tic_poststim);
              Screen('Flip', screen(1));
        vbl = Screen('Flip', screen(2));
        Exit_If_Shift(test,nidaq);
    end

    
    stimtype.tempo(cnt).tempo_stim = tempo;
 
    if exist('scriptName','var') && stimtype.cnt==1
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
end

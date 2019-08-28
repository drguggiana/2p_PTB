
function [stimtype, vbl] = PresentPhiMotion( nidaq, screen, vbl, Param, stimtype )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

for i=1:length(screen)
    Screen('BlendFunction', screen(i), GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
end
cnt=stimtype.cnt;
stimseq=stimtype.stimseq';
test=Param.test;
stimtype.tempo(cnt).toc_poststim = 0;
if Param.Dichoptic
else
    screen(2)=screen(1);
end

        for f = stimtype.screensequence(stimseq(cnt),:) % [1,2] or [2,1]
            tic
            if test == 0
                outputSingleScan(nidaq,stimtype.stim_id);
            end

            for i=1:stimtype.frame_stim    %Ale: # frames per bar
                Screen('FillRect',   screen(f), stimtype.seqcolor(stimseq(cnt)) );    %Ale: fill rectangle with colorBar (black 0 or white 255). Then the mask masktex is applied.
%                 Screen('FillRect',   win, 250 );    %Ale: fill rectangle with colorBar (black 0 or white 255). Then the mask masktex is applied.
                Screen('DrawTexture',screen(f), stimtype.masktex(stimtype.seqbars(stimseq(cnt),f),f) ,[],Param.screenRect(1,:) );  %
                if Param.photodiode
                     Screen('FillRect', screen(1), stimtype.seqcolor(stimseq(cnt)), [0 Param.screenRect(4)-Param.photodiode_size 0+Param.photodiode_size Param.screenRect(4)] );
                end
                Screen('Flip',screen(f));
                Exit_If_Shift(test,nidaq);
            end
%             Screen('FillRect', win, 127);
%             Screen('Flip', win);
%             stimtype.tempo(cnt).toc_bars(f)=toc;
            tempo(1) = toc;
            
             % interpatch
            if f==stimtype.screensequence(1) && stimtype.frame_interpatch ~= 0
                tic
%                 if test == 0
%                     outputSingleScan(nidaq,0);
%                 end
                for ff=1:stimtype.frame_interpatch    %Ale: # frames of blank screen between first and second flashed bar
%                     for i=1:length(screen)
                        Screen('FillRect', screen(f), stimtype.BackgroundLuminance );
                        Screen('Flip', screen(f));
%                     end
                end
                tempo(2) = toc;
                stimtype.tempo(cnt).toc_interpatch=tempo(2);
            end
%             stimtype.tempo(cnt).toc_interpatch=toc;
            Screen('FillRect', screen(f), stimtype.BackgroundLuminance);
            Screen('Flip', screen(f));
        end
        
        %%% post-stim blank screen
        tic;
        if test == 0
           outputSingleScan(nidaq,0);
        end
        for i=1:length(screen)
            Screen('FillRect', screen(i),stimtype.BackgroundLuminance);
        end
        while stimtype.tempo(cnt).toc_poststim < stimtype.poststim_time
            stimtype.tempo(cnt).toc_poststim = toc;
            for i=1:length(screen)
                vbl = Screen('Flip', screen(i));
            end
            Exit_If_Shift(test,nidaq);
        end
%         tempo(3) = stimtype.tempo(cnt).toc_poststim
    
        
    if exist('scriptName','var') && stimtype.cnt==1
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
end

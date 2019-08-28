
function [stimtype, vbl] = PresentPhiMotionBars( nidaq, screen, vbl, Param, stimtype )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));


Screen('BlendFunction', screen, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

Screen('SelectStereoDrawBuffer', screen, 0);
Screen('FillRect', screen, stimtype.BackgroundLuminance);
Screen('SelectStereoDrawBuffer', screen, 1);
Screen('FillRect', screen, stimtype.BackgroundLuminance);
vbl = Screen('Flip', screen);

cnt=stimtype.cnt;
stimseq=stimtype.stimseq';
test=Param.test;
stimtype.tempo(cnt).toc_poststim = 0;


        for s = 1:2
            f = stimtype.screensequence(stimseq(cnt),s);
%         for f = stimtype.screensequence(stimseq(cnt),:) % [1,2] or [2,1]
            if f==0; continue; end;
            %stereo screen where bar is displyed:
            winStereo = f-1;
            % the other stereo screen, which needs to be BackgroundColor:
            if f==1;
                winStereo  = 0;
                winStereo2 = 1;
            elseif f==2
                winStereo  = 1;
                winStereo2 = 0;
            end

            tic

            for i = 1:stimtype.frame_stim    %Ale: # frames per bar
                
                if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
                    Screen('SelectStereoDrawBuffer', screen, winStereo);
                end
                Screen('FillRect',   screen, stimtype.seqcolor(stimseq(cnt)) );    %Ale: fill rectangle with colorBar (black 0 or white 255). Then the mask masktex is applied.
%                 Screen('FillRect',   win, 250 );    %Ale: fill rectangle with colorBar (black 0 or white 255). Then the mask masktex is applied.
                Screen('DrawTexture',screen, stimtype.masktex(stimtype.seqbars(stimseq(cnt),s),f) ,[],Param.screenRect(1,:) );  %
                
                if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
                    Screen('SelectStereoDrawBuffer', screen, winStereo2);
                end
                Screen('FillRect', screen, stimtype.BackgroundLuminance);
                
                if Param.usePhotodiode
                    if ismember(Param.stereoMode,[1 4])
                        Screen('SelectStereoDrawBuffer', screen, Param.Photodiode.screen); 
                    end
                    Screen('FillRect', screen , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
                end
                
                if i==1
                    StimSignal( Param, test, nidaq, stimtype.stim_id )
                end
                
                Screen('DrawingFinished', screen);
                Screen('Flip',screen, vbl + 0.5 * Param.ifi);
                Exit_If_Shift(test,nidaq);
            end
%             Screen('FillRect', win, 127);
%             Screen('Flip', win);
%             stimtype.tempo(cnt).toc_bars(f)=toc;
            tempo(1) = toc;
            
             % interpatch
            if f==stimtype.screensequence(1) && stimtype.frame_interpatch ~= 0
                tic
%                 StimSignal( Param, test, nidaq, 0 )
                for ff=1:stimtype.frame_interpatch    %Ale: # frames of blank screen between first and second flashed bar
%                     for i=1:length(screen)
                        Screen('FillRect', screen, stimtype.BackgroundLuminance );
						if Param.usePhotodiode
                            if ismember(Param.stereoMode,[1 4])
                                Screen('SelectStereoDrawBuffer', screen, Param.Photodiode.screen); 
                            end
							Screen('FillRect', screen, Param.Photodiode.ColorOn, Param.Photodiode.Coord );
                        end
                        Screen('DrawingFinished', screen);
                        Screen('Flip', screen, vbl + 0.5 * Param.ifi);
%                     end
                end
                tempo(2) = toc;
                stimtype.tempo(cnt).toc_interpatch=tempo(2);
            end
%             stimtype.tempo(cnt).toc_interpatch=toc;
            Screen('FillRect', screen, stimtype.BackgroundLuminance);
            Screen('Flip', screen, vbl + 0.5 * Param.ifi);
        end
        
        %%% post-stim blank screen
        tic;
        StimSignal( Param, test, nidaq, 0 )

        Screen('FillRect', screen, stimtype.BackgroundLuminance);

        while stimtype.tempo(cnt).toc_poststim < stimtype.poststim_time
            stimtype.tempo(cnt).toc_poststim = toc;
            vbl = Screen('Flip', screen);
            Exit_If_Shift(test,nidaq);
        end
%         tempo(3) = stimtype.tempo(cnt).toc_poststim
    
        
    if exist('scriptName','var') && stimtype.cnt==1
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
end


function [stimtype, vbl] = PresentRFMapping( nidaq, screen, vbl, Param, stimtype , tar_eye)

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

win = screen(tar_eye);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
% Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


if Param.Dichoptic %&& ismember(Param.stereoMode,[1 4])
%     if     strcmp(inputname(5),'RFM')
%         winStereo = 0;
%         winStereo_notPresented = 1;
%         eye_str = 'left';
%     elseif strcmp(inputname(5),'RFM2')
%         winStereo = 1;
%         winStereo_notPresented = 0;
%         eye_str = 'right';
%     end

%     Screen('SelectStereoDrawBuffer', win, 0);
%     Screen('FillRect', win, stimtype.BackgroundLuminance);
%     Screen('SelectStereoDrawBuffer', win, 1);
%     Screen('FillRect', win, stimtype.BackgroundLuminance);
%     vbl = Screen('Flip', win);
    winStereo = 0; % for Photodiode
    Screen('FillRect', screen(1), stimtype.BackgroundLuminance);
    Screen('FillRect', screen(2), stimtype.BackgroundLuminance);
    vbl = Screen('Flip', screen(1));
    vbl = Screen('Flip', screen(2));
    eye_str = num2str(tar_eye);
elseif Param.Dichoptic==0
    winStereo = 0;
    eye_str = '';
end

cnt = stimtype.cnt;
stimseq = stimtype.stimseq';
test = Param.test;
tempo = [0 0];

disp(['   - Patch nr. ' num2str(stimseq(cnt)) ' ' eye_str ]);


% stimtype.tempo(cnt).toc_poststim = 0;
tic
vblendtime = vbl + stimtype.patch_time - Param.ifi;  %Ale: I put -ifi to obtain an exact stimulus time, I dont know why exactly.
f = 0;

        while vbl < vblendtime
%         for i=1:stimtype.frame_patch    % # frames per bar
            f = f+1;


                if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
                    Screen('SelectStereoDrawBuffer', win, winStereo);
                end
            Screen('FillRect', win, stimtype.colorseq(stimseq(cnt)));    % fill rectangle with black 0 or white 255. Then the mask masktex is applied.
                 
            
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            Screen('DrawTexture', win,stimtype.masktex(stimtype.patchseq(stimseq(cnt))),[],[]);  % patchseq(stimseq(j,s)) returns a random patch (integer from 1 to prod(n_patches)).

%             Disable alpha-blending:
            Screen('Blendfunction', win, GL_ONE, GL_ZERO);
                
%             if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
%                 Screen('SelectStereoDrawBuffer', win, winStereo_notPresented);
%             end
%             Screen('FillRect', win, stimtype.BackgroundLuminance);

                
				if Param.usePhotodiode
                    if ismember(Param.stereoMode,[1 4])
                        Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
                    end
					Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
				end
                Screen('DrawingFinished', win);


                if f==1
                    StimSignal( Param, test, nidaq, stimtype.stim_id )
                end
                % Show it at next retrace:
				vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);
				Exit_If_Shift(test,nidaq);

  

        end
        
        tempo(1) = toc;
        StimSignal( Param, test, nidaq, 0 )
Screen('FillRect', win, stimtype.BackgroundLuminance);
        
        % interpatch
        
        tic 
        if stimtype.frame_interpatch ~= 0
            
%             StimSignal( Param, test, nidaq, 0 )

            vblendtime = vbl + stimtype.interpatch_time - Param.ifi;
            Screen('FillRect', win, stimtype.BackgroundLuminance );
            while vbl < vblendtime
                vbl = Screen('Flip', win);
                Exit_If_Shift(test,nidaq);
            end
        end
        tempo(2) = toc;

        stimtype.tempo{cnt} = tempo;
        
        
    if exist('scriptName','var') && stimtype.cnt==1
        if (length(Param.StimProtocols.RFMapping)>1 && strcmp(inputname(5),'RFM'))...
            || length(Param.StimProtocols.RFMapping)==1
            BackupStimulusScript( scriptPath, scriptName, Param.save_dir )
        end
    end
end
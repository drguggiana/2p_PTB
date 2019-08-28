
function [stimtype, vbl] = PresentRetinotopicPositionStripes( nidaq, win, vbl, Param, stimtype )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
Screen('BlendFunction', win, GL_ONE, GL_ZERO);

if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
    if     strcmp(inputname(5),'RPS')
        winStereo = 0;
    elseif strcmp(inputname(5),'RPS2')
        winStereo = 1;
    end

    Screen('SelectStereoDrawBuffer', win, 0);
    Screen('FillRect', win, stimtype.BackgroundLuminance);
    Screen('SelectStereoDrawBuffer', win, 1);
    Screen('FillRect', win, stimtype.BackgroundLuminance);
    vbl = Screen('Flip', win);
elseif Param.Dichoptic==0
    winStereo = 0;
end

    tic
    cnt = stimtype.cnt;
	test=Param.test;
    if stimtype.directions == 8
        seqangles = [1 4 7 2 5 8 3 6] * 360/stimtype.directions;
    else
        seqangles = randperm(stimtype.directions) * 360/stimtype.directions;
    end
    masktex = stimtype.masktex{stimtype.seqpositions(cnt)};
% 	mask=reshape( stimtype.coord_positions{stimtype.seqpositions(stimtype.cnt)} , 1, 4 );
	tempo=[0 0];
%     tic
% 	vblendtime = vbl + stimtype.stimulus_time - Param.ifi;  %Ale: I put -ifi to obtain an exact stimulus time, I dont know why exactly.
% 	while vbl < vblendtime
        
		for d=1:stimtype.directions
			angle = seqangles(d);
			for f = 1 : fix(stimtype.stimulus_time_fr/stimtype.directions)
				% Update some grating animation parameters:
				% Increment phase by 1 degree:
				stimtype.phase = stimtype.phase + stimtype.phaseincrement;
				% Draw the grating, centered on the screen, with given rotation 'angle',
				% sine grating 'phase' shift and amplitude, rotating via set
				% 'rotateMode'. Note that we pad the last argument with a 4th
				% component, which is 0. This is required, as this argument must be a
				% vector with a number of components that is an integral multiple of 4,
				% i.e. in our case it must have 4 components:
                if Param.Dichoptic && ismember(Param.stereoMode,[1 4])
                    Screen('SelectStereoDrawBuffer', win, winStereo);
                end
				Screen('DrawTexture', win, stimtype.gratingtex, [], [], angle, [], [], [], [], Param.rotateMode, [stimtype.phase, stimtype.freq_cppx, stimtype.amplitude, 0]);
                
                

        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		% Draw mask over grating:
		Screen('DrawTexture', win, masktex, [], []);

    % Disable alpha-blending:
    Screen('Blendfunction', win, GL_ONE, GL_ZERO);
                
				if Param.usePhotodiode
                    if ismember(Param.stereoMode,[1 4])
                        Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
                    end
					Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
				end
                Screen('DrawingFinished', win);
                
                if d==1 && f==1
                    StimSignal( Param, test, nidaq, stimtype.stim_id )
                end
                % Show it at next retrace:
				vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);
				Exit_If_Shift(test,nidaq);
			end
		end
	tempo(1) = toc;
    StimSignal( Param, test, nidaq, 0 )

	tic
	Screen('FillRect', win, stimtype.BackgroundLuminance);
	while  tempo(2) < stimtype.poststim_time;
		tempo(2) = toc;
		Exit_If_Shift(test,nidaq);
        vbl = Screen('Flip', win);  
	end
	stimtype.tempo{stimtype.cnt}= tempo;

    
if exist('scriptName','var') && stimtype.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
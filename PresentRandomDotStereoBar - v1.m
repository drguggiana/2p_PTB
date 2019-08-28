function [RDSB, vbl] = PresentRandomDotStereoBar( nidaq, screen, vbl, Param, RDSB )

% [scriptPath, scriptName] = fileparts(mfilename('fullpath'));


win=screen(1);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
% Screen('SelectStereoDrawBuffer', win, 0);
% Screen('FillRect', win, RDSB.BackgroundLuminance);
% Screen('SelectStereoDrawBuffer', win, 1);
% Screen('FillRect', win, RDSB.BackgroundLuminance);
% vbl = Screen('Flip', win);


test = Param.test;
tempo = [];
cnt = RDSB.cnt;
stimseq        = RDSB.stimseq';
seqIntervals   = RDSB.seqIntervals';
seqangles      = RDSB.seqangles';

disp_ix   = seqIntervals( stimseq(cnt) );
disparity = RDSB.disparity_values(disp_ix);
angle     = seqangles( stimseq(cnt) );

nTexBars = RDSB.nTexBars;

disp([' --> Angle (cart) = ' num2str( mod(angle,360)+90 ) ]);
disp(['    --> Disparity = ' num2str( disparity ) ]);

%% draw and start stimulus
ticstart = tic;
vblendtime = vbl + RDSB.stimulus_time - Param.ifi;
frameNr = 0;
patternNr = 0;

for t = 1 : nTexBars
    
    frameNr = frameNr + 1;
    
    vbl_patternend = vbl + RDSB.pattern_time - Param.ifi;
    if vbl_patternend > vblendtime
        vbl_patternend = vblendtime;
    end
    go_pattern = 1;
    patternNr = patternNr + 1;
    
    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);
	Screen('DrawTexture', win, RDSB.Textures(t,disp_ix,1), [], [], angle, [], [], [], [], Param.rotateMode, []);
	if Param.usePhotodiode
        Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
		Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
    end
    
    % Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 1);
	Screen('DrawTexture', win, RDSB.Textures(t,disp_ix,2), [], [], angle, [], [], [], [], Param.rotateMode, []);

    Screen('DrawingFinished', win);
    
    dontclear = 1;
    while go_pattern == 1
        if frameNr == 1
            StimSignal( Param, test, nidaq, RDSB.stim_id );
        end
        vbl = Screen('Flip', win, vbl + 0.5*Param.ifi, dontclear);
        if vbl > vbl_patternend
            dontclear = 0;
            vbl = Screen('Flip', win, vbl + 0.5*Param.ifi, dontclear);
            go_pattern = 0;
            tempo(patternNr) = toc(ticstart); %#ok<AGROW>
        end
        Exit_If_Shift(test,nidaq);
    end
    
end

StimSignal( Param, test, nidaq, 0 );
% end stimulus

%% post-stimulus
ticpoststim = tic;
Screen('FillRect', win, RDSB.BackgroundLuminance);

for i = 1 : RDSB.frames_poststim
    vbl = Screen('Flip', win);
    Exit_If_Shift(test,nidaq);
end
tempo(end+1) = toc(ticpoststim);
% end post-stimulus
RDSB.tempo{cnt}=tempo;


if exist('scriptName','var') && RDSB.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end

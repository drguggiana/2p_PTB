function [RDS, vbl] = PresentNaturalImageStereogram( nidaq, screen, vbl, Param, RDS, img ) %#ok<INUSL>

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));


win=screen(1);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
Screen('SelectStereoDrawBuffer', win, 0);
Screen('FillRect', win, RDS.BackgroundLuminance);
Screen('SelectStereoDrawBuffer', win, 1);
Screen('FillRect', win, RDS.BackgroundLuminance);
vbl = Screen('Flip', win);

test=Param.test;
tempo=[];
cnt = RDS.cnt;
stimseq      = RDS.stimseq';
seqangles    = RDS.seqangles';
seqIntervals = RDS.seqIntervals';
angle = seqangles(stimseq(cnt));

disp([' --> Angle (cart) = ' num2str(mod(angle,180)) ]);

list_intervals = 1 : RDS.n_disparity_intervals;
whichInterval = list_intervals(seqIntervals(stimseq(cnt)));

bin_ranges_deg = RDS.bin_ranges_deg(whichInterval,:);

seqImages = randperm(1000,RDS.nPatterns);

%% generate and start stimulus
ticstart = tic;
vblendtime = vbl + RDS.stimulus_time - Param.ifi;
frameNr = 0;
patternNr = 0;



while vbl < vblendtime

    frameNr = frameNr + 1;
    
    vbl_patternend = vbl + RDS.pattern_time - Param.ifi;
    if vbl_patternend > vblendtime
        vbl_patternend = vblendtime;
    end
    go_pattern = 1;
    patternNr = patternNr + 1;
    
    IL = img{seqImages(patternNr)};
    [sizeX, sizeY] = size(IL);
    
    disparity_deg = (bin_ranges_deg(2)-bin_ranges_deg(1))*rand(1) + bin_ranges_deg(1);
    dx =  disparity_deg*sin(angle*pi/180) * Param.pixperdeg(1) * sizeX / Param.screenRes(1);
    dy = -disparity_deg*cos(angle*pi/180) * Param.pixperdeg(1) * sizeY / Param.screenRes(2);
    dx = round(dx); % dx>0, displayed image moves rightward
    dy = round(dy); % dy>0, displayed image moves towards bottom

    
    IR = circshift(IL,round(+[dx,dy]/2));
    IL = circshift(IL,round(-[dx,dy]/2));
    
    dstRect = Param.screenRect;
    
    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);
    
    imageTextureL = Screen('MakeTexture', win, IL');

    Screen('DrawTexture', win, imageTextureL, [], dstRect, 0);    
    
    
    % Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 1);
    
    imageTextureR = Screen('MakeTexture', win, IR');

    Screen('DrawTexture', win, imageTextureR, [], dstRect, 0);    

    
    disp([ ' ----> Range disparity [' num2str(bin_ranges_deg) '] - value: ' num2str(disparity_deg,'%2.1f')]);
    
    
    dontclear = 1;
    while go_pattern == 1
        if frameNr == 1
            StimSignal( Param, test, nidaq, RDS.stim_id );
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
Screen('FillRect', win, RDS.BackgroundLuminance);

for i = 1 : RDS.frames_poststim
    vbl = Screen('Flip', win);
    Exit_If_Shift(test,nidaq);
end
tempo(end+1) = toc(ticpoststim);
% end post-stimulus
RDS.tempo{cnt}=tempo;

RDS.seqImages{cnt} = seqImages(1:patternNr);


if exist('scriptName','var') && RDS.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
function [RDSB, vbl] = PresentRandomDotStereoBar( nidaq, screen, vbl, Param, RDSB )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

    test_RDSB_pattern = RDSB.test_RDSB_pattern;

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
stimseq       = RDSB.stimseq';
seqdirections = RDSB.seqdirections';
seqIntervals  = RDSB.seqIntervals';
seqangles     = RDSB.seqangles';

angle_ptb     = seqangles(stimseq(cnt));
angle_cart    = mod(angle_ptb+90,360);
angle_ix      = seqdirections( stimseq(cnt) );
disp_ix       = seqIntervals( stimseq(cnt) );
disparity_deg = RDSB.disparity_values(disp_ix);
disparity_px  = +round( RDSB.disparity_values_px(disp_ix));

nTexBars      = RDSB.nTexBars;

% coord_x_left  = RDSB.coord_x_left (:,:,disp_ix);  % (t,[x1 x2]) 
% coord_x_right = RDSB.coord_x_right(:,:,disp_ix);
% coord_y_left  = RDSB.coord_y_left (:,:,disp_ix);
% coord_y_right = RDSB.coord_y_right(:,:,disp_ix);

mask_bar = RDSB.mask_bar(:,:,:,:,angle_ix,disp_ix); % RDSB.mask(:,:,x,2,a,d) -> mask(:,:,x,2)


%% draw and start stimulus

disp([' --> Angle (cart) = ' num2str( angle_cart ) ]);
disp(['    --> Disparity = ' num2str( disparity_deg ) ]);

ticstart = tic;
vblendtime = vbl + RDSB.stimulus_time - Param.ifi;
frameNr = 0;
patternNr = 0;

for c = 1 : RDSB.n_cycles
for t = 1 : nTexBars
    
    frameNr = frameNr + 1;
    if RDSB.test_RDSB_ComputationTime; tic; end;

    
    vbl_patternend = vbl + RDSB.pattern_time - Param.ifi;
    if vbl_patternend > vblendtime
        vbl_patternend = vblendtime;
    end
    go_pattern = 1;
    patternNr = patternNr + 1;
    
    
    dots.x = (rand(1,RDSB.nDots)-.5)*RDSB.apertureSize(1) + RDSB.center(1);
    dots.y = (rand(1,RDSB.nDots)-.5)*RDSB.apertureSize(2) + RDSB.center(2);
    %convert from degrees to screen pixels
    pixposL.x = dots.x * Param.pixperdeg(1);
    pixposL.y = dots.y * Param.pixperdeg(1); 
    % This generates pixel positions, but they're centered at [0,0]. The last 
    % step for this conversion is to add in the offset for the center of the 
    % screen:
    pixposL.x    = round( pixposL.x + Param.screenRect(3)/2 );
    pixposL.y    = round( pixposL.y + Param.screenRect(4)/2 );
    pixposL.inds = sub2ind([Param.screenRes(2),Param.screenRes(1)], pixposL.y, pixposL.x);

    
    dotsR.x = (rand(1,RDSB.nDots)-.5)*RDSB.apertureSize(1) + RDSB.center(1);
    dotsR.y = (rand(1,RDSB.nDots)-.5)*RDSB.apertureSize(2) + RDSB.center(2);
    
    %convert from degrees to screen pixels
    pixposR.x = dotsR.x * Param.pixperdeg(2);
    pixposR.y = dotsR.y * Param.pixperdeg(2); 
    % This generates pixel positions, but they're centered at [0,0]. The last 
    % step for this conversion is to add in the offset for the center of the 
    % screen:
    pixposR.x    = round( pixposR.x + Param.screenRect(3)/2 );
    pixposR.y    = round( pixposR.y + Param.screenRect(4)/2 );
    pixposR.inds = sub2ind([Param.screenRes(2),Param.screenRes(1)], pixposR.y, pixposR.x);

    
    
    mask = false(Param.screenRes(2), Param.screenRes(1));
    mask(pixposL.inds) = true;
    mask_le = mask .* mask_bar(:,:,t,1);
%     list_pixposL_toCopy = find(mask_le);
    [y_toCopy, x_toCopy] = find(mask_le);
    %
    mask = false(Param.screenRes(2), Param.screenRes(1));
    mask(pixposR.inds) = true;
    mask_re = mask .* mask_bar(:,:,t,2);
    list_pixposR_toDelete = find(mask_re);
%     [y_toDel, x_toDel] = find(mask_re);
    
    list_ix_pixposR_toDelete = find( ismember(pixposR.inds, list_pixposR_toDelete));
    pixposR.x( list_ix_pixposR_toDelete ) = [];
    pixposR.y( list_ix_pixposR_toDelete ) = [];
    
    
    if test_RDSB_pattern == 0
        % append the pixpos copied from left screen to the right screen:
        pixposR.x = [pixposR.x, x_toCopy' + disparity_px];
        pixposR.y = [pixposR.y, y_toCopy' ];
        
    else % test_RDSB_pattern
        % Only the dots within the bar appear, rest is empty:
        pixposR.x = x_toCopy' + disparity_px;
        pixposR.y = y_toCopy' ;
        pixposL.x = x_toCopy';
        pixposL.y = y_toCopy';
    end   
    
%     list_pixposR_toSubstitute = ... % they must be eliminated from pixposR
%         find( pixposR.x >= coord_x_right(t,1) & pixposR.x <= coord_x_right(t,2) &...
%               pixposR.y >= coord_y_right(t,1) & pixposR.y <= coord_y_right(t,2) );
%     
%     list_pixposL_toCopy = ...
%         find( pixposL.x >= coord_x_left(t,1)  & pixposL.x <= coord_x_left(t,2)  &...
%               pixposL.y >= coord_y_left(t,1)  & pixposL.y <= coord_y_left(t,2)  );
% 
%     pixposR.x(list_pixposR_toSubstitute) = [];
%     pixposR.y(list_pixposR_toSubstitute) = [];
%     if test_RDSB_pattern == 0
%         % append the pixpos copied from left screen to the right screen:
%         pixposR.x = [pixposR.x, pixposL.x(list_pixposL_toCopy) - disparity_px];
%         pixposR.y = [pixposR.y, pixposL.y(list_pixposL_toCopy)];
%     else % test_RDSB_pattern
%         % Only the dots within the bar appear, rest is empty:
%         pixposR.x = pixposL.x(list_pixposL_toCopy) - disparity_px;
%         pixposR.y = pixposL.y(list_pixposL_toCopy);
%         pixposL.x = pixposL.x(list_pixposL_toCopy);
%         pixposL.y = pixposL.y(list_pixposL_toCopy);
%     end
    
    
    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);
    
    Screen('DrawDots',win, [pixposL.x; pixposL.y], RDSB.size_px(1), RDSB.color{1},[0,0],1);
    
    if Param.usePhotodiode
        Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
		Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
    end
    
    
    % Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 1);    
    
    Screen('DrawDots',win, [pixposR.x; pixposR.y], RDSB.size_px(2), RDSB.color{1},[0,0],1);


    Screen('DrawingFinished', win);
    
    
    if RDSB.test_RDSB_ComputationTime; toc; end;

    
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

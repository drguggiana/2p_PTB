function [RDSE, vbl] = PresentRandomDotStereoEdge( nidaq, screen, vbl, Param, RDSE )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));


win=screen(1);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
% Screen('SelectStereoDrawBuffer', win, 0);
% Screen('FillRect', win, RDSE.BackgroundLuminance);
% Screen('SelectStereoDrawBuffer', win, 1);
% Screen('FillRect', win, RDSE.BackgroundLuminance);
% vbl = Screen('Flip', win);

test  = Param.test;
tempo = [];
cnt   = RDSE.cnt;
stimseq       = RDSE.stimseq';
seqdirections = RDSE.seqdirections';
seqIntervals  = RDSE.seqIntervals';
seqangles     = RDSE.seqangles';

angle_ptb     = seqangles(stimseq(cnt));
angle_cart    = mod(angle_ptb+90,360);
angle_ix      = seqdirections( stimseq(cnt) );
disp_ix       = seqIntervals( stimseq(cnt) );
disparity_deg = RDSE.disparity_values(disp_ix);
if angle_cart > 180
    disphalf_px  = -round( RDSE.disparity_values_px(disp_ix)/2 );
else
    disphalf_px  = +round( RDSE.disparity_values_px(disp_ix)/2);
end

nTexBars = RDSE.nTexBars;

mask_edge = RDSE.mask_edge(:,:,:,angle_ix);
mask_gap  = RDSE.mask_gap (:,:,:,angle_ix,disp_ix);


%% generate and start stimulus

disp([' --> Angle (cart) = ' num2str( angle_cart ) ]);
disp(['    --> Disparity = ' num2str( disparity_deg ) ]);

ticstart = tic;
vblendtime = vbl + RDSE.stimulus_time - Param.ifi;
frameNr = 0;
patternNr = 0;

for c = 1 : RDSE.n_cycles
for t = 1 : nTexBars
    
    if RDSE.test_RDSE_ComputationTime; tic; end;
    frameNr = frameNr + 1;
    
    vbl_patternend = vbl + RDSE.pattern_time - Param.ifi;
    if vbl_patternend > vblendtime
        vbl_patternend = vblendtime;
    end
    go_pattern = 1;
    patternNr = patternNr + 1;
    
    
    %%% Dots at left screen
    
    dots.x = (rand(1,RDSE.nDots)-.5)*RDSE.apertureSize(1) + RDSE.center(1);
    dots.y = (rand(1,RDSE.nDots)-.5)*RDSE.apertureSize(2) + RDSE.center(2);
    %convert from degrees to screen pixels
    pixposL.x =  dots.x * Param.pixperdeg(1) ;
    pixposL.y =  dots.y * Param.pixperdeg(1) ; 
    % This generates pixel positions, but they're centered at [0,0]. The last 
    % step for this conversion is to add in the offset for the center of the 
    % screen:
    % It's important to round here! Otherwise you won't find the dots to
    % copy:
    pixposL.x    = round( pixposL.x + Param.screenRect(3)/2 );
    pixposL.y    = round( pixposL.y + Param.screenRect(4)/2 );
    pixposL.inds = sub2ind([Param.screenRes(2),Param.screenRes(1)], pixposL.y, pixposL.x);
    % the following approach to get pixposL is NOT faster:
    % %     a = 1; b = prod(Param.screenRes);
    % %     pixposL.inds = round( (b-a).*rand(RDSE.nDots,1) + a);
    % %     [pixposL.y,pixposL.x] = ind2sub([Param.screenRes(2),Param.screenRes(1)],pixposL.inds);
    
    %%% Dots at right screen
        
    mask = false(Param.screenRes(2), Param.screenRes(1));
    mask(pixposL.inds) = true;
    mask_le = mask .* mask_edge(:,:,t);
    [pixposR.y, pixposR.x] = find(mask_le);
    pixposR.x = pixposR.x - disphalf_px;
    %
    mask_re = mask .* ~mask_edge(:,:,t) ;
    [y_re, x_re] = find(mask_re);
    pixposR.x = [pixposR.x; x_re + disphalf_px]';
    pixposR.y = [pixposR.y; y_re ]';
    
    % the following approach to get pixposR is slower, it's better to use
    % binary masks than indices:
%     edge_le = find( mask_edge(:,:,t));
%     edge_re = find(~mask_edge(:,:,t));
%     tic
%     inds_le = intersect( edge_le,pixposL.inds);
%     [pixposR.y,pixposR.x] = ind2sub([Param.screenRes(2),Param.screenRes(1)], inds_le);
%     pixposR.x = pixposR.x - disparity_px/2;
%     inds_re = intersect( edge_re,pixposL.inds);
%     [y,x] = ind2sub([Param.screenRes(2),Param.screenRes(1)], inds_re);
%     pixposR.x = [pixposR.x; x + disparity_px/2]';
%     pixposR.y = [pixposR.y; y ]';
%     toc
    
     
%     list_le = find(ismember(pixposL.x, x_le{t}) & ismember(pixposL.y, y_le{t}));
%     list_re = find(ismember(pixposL.x, x_re{t}) & ismember(pixposL.y, y_re{t}));
%     list_le = find( ismember([pixposL.x',pixposL.y'], xy_le{t}, 'rows') );
%     list_re = find( ismember([pixposL.x',pixposL.y'], xy_re{t}, 'rows') );
%     list_le = intersect(pixposL.inds, ind_le{t});
%     list_re = intersect(pixposL.inds, ind_re{t});
% 
%     pixposR.x = pixposL.x(list_le) - disparity_px/2;
%     pixposR.y = pixposL.y(list_le);
%     %
%     pixposR.x = [pixposR.x, pixposL.x(list_re) + disparity_px/2];
%     pixposR.y = [pixposR.y, pixposL.y(list_re)];
    

    
    if ismember(angle_cart, [0, 180])
        % With horizontal edge there is no gap at the level of the edge, but we
        % create gaps at the margins of the screen, so we circularly move back the dots
        % that fall outside:
        % /!\ This correction is implemented only for horizontal edge for now;
        % oblique edge is bit more complicated, not crucial for now...
        pixposR.x( pixposR.x > Param.screenRect(3)) = pixposR.x( pixposR.x > Param.screenRect(3) ) - Param.screenRect(3);
        pixposR.x( pixposR.x < 1                  ) = pixposR.x( pixposR.x < 1 ) + Param.screenRect(3);
    else
        
        % now remains the gap
        %draw random dots again:
        dots_gap.x = (rand(1,RDSE.nDots)-.5)*RDSE.apertureSize(1) + RDSE.center(1);
        dots_gap.y = (rand(1,RDSE.nDots)-.5)*RDSE.apertureSize(2) + RDSE.center(2);
        %convert from degrees to screen pixels
        pixpos_gap.x =  dots_gap.x * Param.pixperdeg(1) ;
        pixpos_gap.y =  dots_gap.y * Param.pixperdeg(1) ; 

        pixpos_gap.x = round( pixpos_gap.x + Param.screenRect(3)/2 );
        pixpos_gap.y = round( pixpos_gap.y + Param.screenRect(4)/2 );
        pixpos_gap.inds = sub2ind([Param.screenRes(2),Param.screenRes(1)], pixpos_gap.y, pixpos_gap.x);

        %         list_gap = find(ismember(pixpos_gap.x, x_gap{t}) & ismember(pixpos_gap.y, y_gap{t}))';
%         [~,list_gap] = ismember([pixpos_gap.x',pixpos_gap.y'], ind_gap{t}, 'rows');
        
        mask = false(Param.screenRes(2), Param.screenRes(1));
        mask(pixpos_gap.inds) = true;
        mask2 = mask .* mask_gap(:,:,t);
        [y_g, x_g] = find(mask2);
        pixposR.x = [pixposR.x, x_g' ];
        pixposR.y = [pixposR.y, y_g' ];
        
%         pixposR.x = [pixposR.x, pixpos_gap.x(list_gap)];
%         pixposR.y = [pixposR.y, pixpos_gap.y(list_gap)];
        
%     pixposR.x =  pixpos_gap.x(list_gap);
%     pixposR.y =  pixpos_gap.y(list_gap);

    end
    

    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);
    
    Screen('DrawDots',win, [pixposL.x; pixposL.y], RDSE.size_px(1), RDSE.color{1},[0,0],1);
    
    % Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 1);
    
    Screen('DrawDots',win, [pixposR.x; pixposR.y], RDSE.size_px(2), RDSE.color{1},[0,0],1);
    
    Screen('DrawingFinished', win);
    
    if RDSE.test_RDSE_ComputationTime; toc; end;
    
    dontclear = 1;
    while go_pattern == 1
        if frameNr == 1
            StimSignal( Param, test, nidaq, RDSE.stim_id );
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
if RDSE.poststim_time > 0
    Screen('FillRect', win, RDSE.BackgroundLuminance);

    for i = 1 : RDSE.frames_poststim
        vbl = Screen('Flip', win);
        Exit_If_Shift(test,nidaq);
    end
end
tempo(end+1) = toc(ticpoststim);
% end post-stimulus
RDSE.tempo{cnt}=tempo;


if exist('scriptName','var') && RDSE.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
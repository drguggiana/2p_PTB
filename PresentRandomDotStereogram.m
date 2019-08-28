function [RDS, vbl] = PresentRandomDotStereogram( nidaq, screen, vbl, Param, RDS ) %#ok<INUSL>

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

% Now we'll define a random position within the aperture for each of the
% DM. 'DM.x' and 'DM.y' will hold the x and y positions for each dot.
% What's the logic here? 'rand(DM.nDots,1)' generates a column vector
% DM.nDots long of random numbers between 0 and 1. To change that range
% to fit the aperture we subtract .5 from those numbers, multiply them by
% the aperture size, and add the center offset. Get it?
% dots.x = (rand(1,DM.nDots)-.5)*DM.apertureSize(1) + DM.center(1);
% dots.y = (rand(1,DM.nDots)-.5)*DM.apertureSize(2) + DM.center(2);

% Each dot will have a integer value 'life' which is how many frames the %
% dot has been going. The starting 'life' of each dot will be a random %
% number between 0 and DM.lifetime-1 so that they don't all 'die' on the
% same frame:
% dots.life = ceil(rand(1,DM.nDots)*DM.lifetime_fr);
dots.life = zeros(1,RDS.nDots);

% Animation is performed by updating the dot position on each frame and
% re-drawing the frame. We need to know the frame-rate of our monitor so
% that we can calculate how much we need to change the dot positions on
% each frame. Fortuantely, our 'OpenWindow' function appends the field
% 'frameRate' to the 'display' structure.
% The distance traveled by a dot (in degrees) is the speed (degrees/second)
% divided by the frame rate (frames/second). The units cancel, leaving
% degrees/frame which makes sense. Basic trigonometry (sines and cosines)
% allows us to determine how much the changes in the x and y position. So
% the x and y position changes, which we'll call dx and dy, can be
% calculated by:
% phi_jump = (DM.phi_jump(2)-DM.phi_jump(1))*rand(1) + DM.phi_jump(1);
% dx =  phi_jump*sin(angle*pi/180);
% dy = -phi_jump*cos(angle*pi/180);

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
    
    
    dots.x = (rand(1,RDS.nDots)-.5)*RDS.apertureSize(1) + RDS.center(1);
    dots.y = (rand(1,RDS.nDots)-.5)*RDS.apertureSize(2) + RDS.center(2);
    %convert from degrees to screen pixels
    pixpos.x = dots.x * Param.pixperdeg(1);
    pixpos.y = dots.y * Param.pixperdeg(1); 
    % This generates pixel positions, but they're centered at [0,0]. The last 
    % step for this conversion is to add in the offset for the center of the 
    % screen:
    pixpos.x = pixpos.x + Param.screenRect(3)/2;
    pixpos.y = pixpos.y + Param.screenRect(4)/2;

    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);
    
    Screen('DrawDots',win, [pixpos.x;pixpos.y], RDS.size_px(1), RDS.color_list{1},[0,0],1);
    
    
    
    % Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 1);
       
    
    disparity_deg = (bin_ranges_deg(2)-bin_ranges_deg(1))*rand(1) + bin_ranges_deg(1);
    dx =  disparity_deg*sin(angle*pi/180);
    dy = -disparity_deg*cos(angle*pi/180);
    
    disp([ ' ----> Range disparity [' num2str(bin_ranges_deg) '] - value: ' num2str(disparity_deg,'%2.1f')]);
    
    %update the dot position
    dotsR.x = dots.x + dx;
    dotsR.y = dots.y + dy;

    
    %move the dots that are outside the aperture back one aperture
    %width.
    l = RDS.aperture.l;
    r = RDS.aperture.r;
    b = RDS.aperture.b;
    t = RDS.aperture.t;
    dotsR.x(dots.x<l) = dotsR.x(dots.x<l) + RDS.apertureSize(1);
    dotsR.x(dots.x>r) = dotsR.x(dots.x>r) - RDS.apertureSize(1);
    dotsR.y(dots.y<b) = dotsR.y(dots.y<b) + RDS.apertureSize(2);
    dotsR.y(dots.y>t) = dotsR.y(dots.y>t) - RDS.apertureSize(2);
    
%     %increment the 'life' of each dot
%     dots.life = dots.life+1;
% 
%     %find the 'dead' dots
%     deadDots = mod(dots.life,RDS.lifetime_fr)==0;
% 
%     %replace the positions of the dead dots to a random location
%     dots.x(deadDots) = (rand(1,sum(deadDots))-.5)*RDS.apertureSize(1) + RDS.center(1);
%     dots.y(deadDots) = (rand(1,sum(deadDots))-.5)*RDS.apertureSize(2) + RDS.center(2);
    
    %convert from degrees to screen pixels
    pixpos.x = dotsR.x * Param.pixperdeg(2);
    pixpos.y = dotsR.y * Param.pixperdeg(2); 
    % This generates pixel positions, but they're centered at [0,0]. The last 
    % step for this conversion is to add in the offset for the center of the 
    % screen:
    pixpos.x = pixpos.x + Param.screenRect(3)/2;
    pixpos.y = pixpos.y + Param.screenRect(4)/2;

    Screen('DrawDots',win, [pixpos.x;pixpos.y], RDS.size_px(2), RDS.color_list{2},[0,0],1);
    
    
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


if exist('scriptName','var') && RDS.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
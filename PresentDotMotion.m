function [DM, vbl] = PresentDotMotion( nidaq, screen, vbl, Param, DM )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if     Param.Dichoptic && Param.stereoMode==0 
    win=screen(1);
    winStereo = 0; % for Photodiode
    Screen('FillRect', screen(1), DM.BackgroundLuminance);
    Screen('FillRect', screen(2), DM.BackgroundLuminance);
    vbl = Screen('Flip', screen(1));
    vbl = Screen('Flip', screen(2));
elseif Param.Dichoptic && Param.stereoMode==4
    winStereo = 0;
    win=screen(1);
    Screen('SelectStereoDrawBuffer', win, 0);
    Screen('FillRect', win, DM.BackgroundLuminance);
    Screen('SelectStereoDrawBuffer', win, 1);
    Screen('FillRect', win, DM.BackgroundLuminance);
    vbl = Screen('Flip', win);
else
    win=screen(1);
    winStereo = 0;
    vbl = Screen('Flip', win);
end

% Screen('BlendFunction', win, GL_ONE, GL_ZERO);
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
test=Param.test;
tempo=[0,0];
seqangles = DM.seqangles';
angle = seqangles(DM.cnt);

% Now we'll define a random position within the aperture for each of the
% DM. 'DM.x' and 'DM.y' will hold the x and y positions for each dot.
% What's the logic here? 'rand(DM.nDots,1)' generates a column vector
% DM.nDots long of random numbers between 0 and 1. To change that range
% to fit the aperture we subtract .5 from those numbers, multiply them by
% the aperture size, and add the center offset. Get it?
dots.x = (rand(1,DM.nDots)-.5)*DM.apertureSize(1) + DM.center(1);
dots.y = (rand(1,DM.nDots)-.5)*DM.apertureSize(2) + DM.center(2);

% Each dot will have a integer value 'life' which is how many frames the %
% dot has been going. The starting 'life' of each dot will be a random %
% number between 0 and DM.lifetime-1 so that they don't all 'die' on the
% same frame:
% dots.life = ceil(rand(1,DM.nDots)*DM.lifetime_fr);
dots.life = zeros(1,DM.nDots);

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
dx =  DM.speed*sin(angle*pi/180)/(1/Param.ifi);
dy = -DM.speed*cos(angle*pi/180)/(1/Param.ifi);



for i = 1 : DM.frames_stim
    tic

    %convert from degrees to screen pixels
    pixpos.x = dots.x * Param.pixperdeg(1);
    pixpos.y = dots.y * Param.pixperdeg(1); 
    % This generates pixel positions, but they're centered at [0,0]. The last 
    % step for this conversion is to add in the offset for the center of the 
    % screen:
    % 
    pixpos.x = pixpos.x + Param.screenRect(3)/2;
    pixpos.y = pixpos.y + Param.screenRect(4)/2;
    
    if Param.Dichoptic && Param.stereoMode==4
        Screen('SelectStereoDrawBuffer', win, winStereo);
    end
    Screen('DrawDots',win, [pixpos.x;pixpos.y], DM.size_px, DM.color_list, [0,0], 1);
    % 0 (default) squares, 1 circles (with anti-aliasing), 2 circles (with high-quality anti-aliasing, if supported by your hardware)

    %update the dot position
    dots.x = dots.x + dx;
    dots.y = dots.y + dy;

    %move the dots that are outside the aperture back one aperture
    %width.
    l = DM.aperture.l;
    r = DM.aperture.r;
    b = DM.aperture.b;
    t = DM.aperture.t;
    dots.x(dots.x<l) = dots.x(dots.x<l) + DM.apertureSize(1);
    dots.x(dots.x>r) = dots.x(dots.x>r) - DM.apertureSize(1);
    dots.y(dots.y<b) = dots.y(dots.y<b) + DM.apertureSize(2);
    dots.y(dots.y>t) = dots.y(dots.y>t) - DM.apertureSize(2);

    %increment the 'life' of each dot
    dots.life = dots.life+1;

    %find the 'dead' dots
    deadDots = mod(dots.life,DM.lifetime_fr)==0;

    %replace the positions of the dead dots to a random location
    dots.x(deadDots) = (rand(1,sum(deadDots))-.5)*DM.apertureSize(1) + DM.center(1);
    dots.y(deadDots) = (rand(1,sum(deadDots))-.5)*DM.apertureSize(2) + DM.center(2);

    if Param.usePhotodiode
        Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
		Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
    end
    
    Screen('DrawingFinished', win);
    
    if i == 1 
		StimSignal( Param, test, nidaq, DM.stim_id )
    end
    
    Screen('Flip', win, vbl + 0.5 * Param.ifi);
    Exit_If_Shift(test,nidaq);
end
tempo(1) = toc;
% end stimulus

if test == 0
	StimSignal( Param, test, nidaq, 0 )
end

% post-stimulus
Screen('FillRect', win, DM.BackgroundLuminance);

for i = 1:DM.frames_poststim
    vbl = Screen('Flip', win);
    Exit_If_Shift(test,nidaq);
end
tempo(2) = toc;
% end post-stimulus
DM.tempo{DM.cnt}=tempo;

if exist('scriptName','var') && DM.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
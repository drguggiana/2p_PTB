function Photodiode = SetPhotodiode(screenRect)
% Set parameters of photodiode square presentation

% set size of PD square:
Photodiode.Size     =  50; %in pixels

% set color of square during actual stimulus presentation:
Photodiode.ColorOn  = 255; 

% set color of square during interstimulus interval:
Photodiode.ColorOff =   0;

% set in which screen (l/r) present the square, in case of stereoMode==4
Photodiode.screen = 0; % left screen
% Photodiode.screen = 1; % right screen

% set coordinates of PD rectangle:
Photodiode.Coord    = [	0;... left
						screenRect(4)-Photodiode.Size;... top
						0+Photodiode.Size;... right
						screenRect(4) ]; % bottom

% Clear the workspace and the screen
sca;
close all;
clearvars;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = 1;
black = 0;
grey = white / 2;

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Here we calculate the radial distance from the center of the screen to
% the X and Y edges
xRadius = windowRect(3) / 2;
yRadius = windowRect(4) / 2;

% Screen resolution in Y
screenYpix = windowRect(4);

% Number of white/black circle pairs
rcycles = 5;

% Number of white/black angular segment pairs (integer)
tcycles = 12;

% Now we make our checkerboard pattern
xylim = 2 * pi * rcycles;
[x, y] = meshgrid(-xylim: 2 * xylim / (screenYpix - 1): xylim,...
    -xylim: 2 * xylim / (screenYpix - 1): xylim);
at = atan2(y, x);
checks = ((1 + sign(sin(at * tcycles) + eps)...
    .* sign(sin(sqrt(x.^2 + y.^2)))) / 2) * (white - black) + black;
circle = x.^2 + y.^2 <= xylim^2;
checks = circle .* checks + grey * ~circle;

% Now we make this into a PTB texture
radialCheckerboardTexture  = Screen('MakeTexture', window, checks);

% Draw our texture to the screen
Screen('DrawTexture', window, radialCheckerboardTexture);

% Flip to the screen
Screen('Flip', window);

% Wait for a keypress
KbStrokeWait;

% Clear up and leave the building
sca;

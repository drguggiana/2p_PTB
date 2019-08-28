function [ G, xSteps, ySteps, Specs ] = SimpleGrating( Angle, Options )
% function [ G, xSteps, ySteps, Specs ] = SimpleGrating( Angle, Options )
% 
% Creates a single grating with the properties set in options.
%
% inputs:
% - Angle:      Orientation angle of the grating
% - Options:    Empty by default
%   - .spatialF:    Spatial frequency in cycles per degree (default=0.05)
%   - .temporalF:   Temporal frequency in cycles per second (default=2)
%   - .refreshRate: Refresh rate of the screen (default=60)
%   - .contrastmin: Minimum intensity of the grating (0-255) (default=0)
%   - .contrastmax: Maximum intensity of the grating (0-255) (default=255)
%   - .shape:       'Sinusoid' or 'square' (default='square')
%   - .xsize:       Number of pixels of the grating in the X dimension (default=800)
%   - .ysize:       Number of pixels of the grating in the Y dimension (default=600)
%   - .scrwidth:    Width of the screen in cm (default=34)
%   - .scrheight:   Height of the screen in cm (default=27.4)
%   - .scrdist:     Distance of the eye to the screen in cm (default=16)
%   - .scrx:        Number of pixels of the screen in the X dimension (default=640)
%   - .scry:        Number of pixels of the screen in the Y dimension (default=480)
%
% output:
% - G:      (ysize*xsize) matrix containing grating
% - xSteps: steps of grating displacement on x axes
% - ySteps: steps of grating displacement on y axes
% - Specs:  Specifications of the screen
%     .ScrWidthDegree:      Screen width in retinal degrees
%     .ScrHeightDegree:     Screen height in retinal degrees
%     .PixelWidthDegree:    Pixel width in retinal degrees
%     .PixelHeightDegree:   Pixel height in retinal degrees
%
% written by Pieter Goltstein
%


%% set options to defaults if they are not set

if isempty( Options )
    Options = struct;
end
if ~isfield( Options, 'spatialF' )
    Options.spatialF = 0.05; % Cycles per degree
end
if ~isfield( Options, 'temporalF' )
    Options.temporalF = 2; % Cycles per second
end
if ~isfield( Options, 'refreshRate' )
    Options.refreshRate = 60; % Hz
end
if ~isfield( Options, 'contrastmin' )
    Options.contrastmin = 0;
end
if ~isfield( Options, 'contrastmax' )
    Options.contrastmax = 255;
end
if ~isfield( Options, 'shape' )
    Options.shape = 'square';
end
if ~isfield( Options, 'scrwidth' )
    Options.scrwidth = 34; % cm
end
if ~isfield( Options, 'scrheight' )
    Options.scrheight = 27.3; % cm
end
if ~isfield( Options, 'scrdist' )
    Options.scrdist = 16; % cm
end
if ~isfield( Options, 'scrx' )
    Options.scrx = 640; % pixels
end
if ~isfield( Options, 'scry' )
    Options.scry = 480; % pixels
end
if ~isfield( Options, 'xsize' )
    Options.xsize = round(Options.scrx*1.5); % pixels
end
if ~isfield( Options, 'ysize' )
    Options.ysize = round(Options.scry*1.5); % pixels
end


%% calculate additional size information
ScrWidthDegree  = 2 * atand( (0.5*Options.scrwidth) /Options.scrdist);
ScrHeightDegree = 2 * atand( (0.5*Options.scrheight)/Options.scrdist);

PixelWidthDegree  = ScrWidthDegree  / Options.scrx;
PixelHeightDegree = ScrHeightDegree / Options.scry;

xsizeDegree = Options.xsize * PixelWidthDegree;
ysizeDegree = Options.ysize * PixelHeightDegree;

ContrastDiff = (Options.contrastmax - Options.contrastmin);

Specs.ScrWidthDegree = ScrWidthDegree;
Specs.ScrHeightDegree = ScrHeightDegree;
Specs.PixelWidthDegree = PixelWidthDegree;
Specs.PixelHeightDegree = PixelHeightDegree;


%% Calculate grating

% Get angle in radians
AngleRad = Angle*pi/180;

% create meshgrid of the visual field in retinal degree's
[ xDeg, yDeg ] = meshgrid( linspace(-0.5*xsizeDegree, 0.5*xsizeDegree, Options.xsize ), ...
                linspace(-0.5*ysizeDegree, 0.5*ysizeDegree, Options.ysize ) );

% calculate spatial frequency variable
xFreq = cos(AngleRad) * Options.spatialF * 2 * pi;
yFreq = sin(AngleRad) * Options.spatialF * 2 * pi;

% calculate grating
G = Options.contrastmin + ContrastDiff * ( ( 1 + sin( (xFreq * xDeg) + (yFreq*yDeg) ) ) / 2 );

% check if we need a square wave grating or a sin grating
if strcmpi(Options.shape, 'square') == 1
    SQ = ones(size(G)) .* Options.contrastmin;
    SQ(G >= (Options.contrastmin+(0.5*ContrastDiff))) = Options.contrastmax;
    G = SQ;
end

%% set limits to intensity to prevent psychtoolbox from getting upset..
G(G<0) = 0;
G(G>255) = 255;

%% calculate vector for x and y movement steps
xCyclePx = 1/(Options.spatialF * PixelWidthDegree);
yCyclePx = 1/(Options.spatialF * PixelHeightDegree);
xSteps = linspace( 0, xCyclePx, round((Options.refreshRate/Options.temporalF) + 1) );
ySteps = linspace( 0, yCyclePx, round((Options.refreshRate/Options.temporalF) + 1) );
xSteps = xSteps(1:end-1);
ySteps = ySteps(1:end-1);
xSteps = round(cos(AngleRad).*xSteps);
ySteps = round(sin(AngleRad).*ySteps);

% make sure that steps have only positive values
xSteps = xSteps - min(xSteps);
ySteps = ySteps - min(ySteps);





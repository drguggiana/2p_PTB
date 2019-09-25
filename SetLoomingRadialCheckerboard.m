	
% set some Radial Checkerboard parameters, saved in the structure 'LO'.
% This creates a radial checkerboarf that will be set as the background for
% the looming stim
function [ LO, Param ] = SetLoomingRadialCheckerboard( win, LO, Param )
    [scriptPath, scriptName] = fileparts(mfilename('fullpath'));

    pixperdeg = Param.pixperdeg(1);

	LO.gratingRect = [Param.screenRect(1), Param.screenRect(2), Param.screenRect(3), Param.screenRect(4)];
    LO.gratingsize = [Param.screenRect(3), Param.screenRect(4)]; % x,y
    
    %% build a radial checkerboard
    
    % Define black and white
    white = 1;
    black = 0;
    grey = white / 2;
    
    % Here we calculate the radial distance from the center of the screen to
    % the X and Y edges
    xRadius = LO.gratingsize(1) / 2;
    yRadius = LO.gratingsize(2) / 2;

    % Screen resolution in Y
    screenYpix = LO.gratingRect(4);

    % Number of white/black circle pairs
    rcycles = 8;

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
    LO.radialChecker = Screen('MakeTexture', window, checks);
    

    %% Expansion parameters
	%get the number of conditions
    n_speeds = length(LO.expansion_speeds);
    cond_num = n_speeds;
    
    %assemble a matrix with all the parameters to be varied within a trial
    LO.paramorder = zeros(LO.n_reps, cond_num, 2);
    
    for r = 1 : LO.n_reps
        switch Param.seqmode
            case 'random'
                [x,y] = ind2sub(n_speeds,randperm(cond_num));
                LO.paramorder(r,:,:) = permute(cat(1,x,y),[3 2 1]);
            case 'sequential'
                [x,y] = ind2sub(n_speeds,1:cond_num);
                LO.paramorder(r,:,:) = cat(2,x,y);
        end
    end
    
    %define the expansion speeds in terms of pix/s
    LO.seqspeeds = pixperdeg .* LO.expansion_speeds(LO.paramorder(:,:,1)) .* Param.ifi;
    LO.seqtimes = LO.time_perstimulus(LO.paramorder(:,:,1));
    Param.stimSeq(end+1,1) = { repmat(LO.stim_id, 1, n_speeds * LO.n_reps) };
    LO.cnt=1;
    
    % define the 
   
if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
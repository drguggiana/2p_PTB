	
% set some DriftingGratings parameters, saved in the structure 'LO'.
function [ LO, Param ] = SetLooming( win, LO, Param )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

    pixperdeg = Param.pixperdeg(1);
%     LO.freq_cppx = LO.spacfreq/pixperdeg;
   

	% Build a procedural sine grating texture for a grating with a support of
	% gratingsize(1) x gratingsize(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
	% you need to have the mex file CreateProceduralGratingMod and the GSGL
	% shaders.
	LO.gratingRect = [Param.screenRect(1), Param.screenRect(2), Param.screenRect(3), Param.screenRect(4)];
    LO.gratingsize = [Param.screenRect(3), Param.screenRect(4)]; % x,y
    
% 	if LO.sinwave
% 		LO.gratingtex = CreateProceduralGratingMod(win, 'sin', LO.gratingsize(1), LO.gratingsize(2), [LO.amplitude LO.amplitude LO.amplitude 0], [LO.amplitude LO.amplitude LO.amplitude 0]);
%         % default gratingtex, with this you cannot
%         % set square but only sine. On the other hand you don't need extra
%         % mex files and GSGLShaders:
% %         LO.gratingtex = CreateProceduralSineGrating(win, LO.gratingsize(1), LO.gratingsize(2), [0.5 0.5 0.5 0]);
% 	else
% 		LO.gratingtex = CreateProceduralGratingMod(win, 'square', LO.gratingsize(1), LO.gratingsize(2), [LO.amplitude LO.amplitude LO.amplitude LO.amplitude], [LO.amplitude LO.amplitude LO.amplitude 0]);
% 	end
	
    %% build a radial checkerboard
    if LO.useChecker
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
        LO.radialChecker = Screen('MakeTexture', win, checks);
    end

	%get the number of conditions
    n_speeds = length(LO.expansion_speeds);
    n_colors = size(LO.colors,1);
    
    %number of conditions
    cond_num = n_speeds*n_colors;
    %assemble a matrix with all the parameters to be varied within a trial
    LO.paramorder = zeros(LO.n_reps,cond_num,2);
    
    for r = 1 : LO.n_reps
        switch Param.seqmode
            case 'random'
                [x,y] = ind2sub([n_speeds,n_colors],randperm(cond_num));
                LO.paramorder(r,:,:) = permute(cat(1,x,y),[3 2 1]);
            case 'sequential'
                [x,y] = ind2sub([n_speeds,n_colors],1:cond_num);
                LO.paramorder(r,:,:) = cat(2,x,y);
        end
    end
    
    %define the expansion speeds in terms of pix/s
    LO.seqspeeds = pixperdeg.*LO.expansion_speeds(LO.paramorder(:,:,1)).*Param.ifi;
    LO.seqcolors = LO.colors(LO.paramorder(:,:,2)',:);
    LO.seqtimes = LO.time_perstimulus(LO.paramorder(:,:,1));
    Param.stimSeq(end+1,1) = { repmat(LO.stim_id, 1,n_colors*n_speeds*LO.n_reps) };
    LO.cnt=1;
    
%     for r = 1 : LO.n_reps
%         switch Param.seqmode
%             case 'random'
%                 
%                 LO.seqspeedorder(r,:) = randperm(LO.n_speeds);
%             case 'sequential'
%                 LO.seqspeedorder(r,:) = [1:LO.n_speeds];
%         end
%     end
%     
%     %define the expansion speeds in terms of pix/s
%     LO.seqspeeds = pixperdeg.*LO.expansion_speeds(LO.seqspeedorder).*Param.ifi;
%     LO.seqtimes = LO.time_perstimulus(LO.seqspeedorder);
%     Param.stimSeq(end+1,1) = { repmat(LO.stim_id, 1,LO.n_speeds*LO.n_reps) };
%     LO.cnt=1;


    
if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
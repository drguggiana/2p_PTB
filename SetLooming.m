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
        white = WhiteIndex(Param.tar_eye);
        black = BlackIndex(Param.tar_eye);
        grey = white / 2;

        % Screen resolution in Y
        screenYpix = Param.screenRect(4);

        % Number of white/black circle pairs
        LO.degpercycle = 5;
        LO.rcycles = floor(LO.screen_angularSize/LO.degpercycle);

        % Number of white/black angular segment pairs (integer)
        LO.tcycles = 12;

        % Now we make our checkerboard pattern
        xylim = 2 * pi * LO.rcycles;
        [x, y] = meshgrid(-xylim: 2 * xylim / (screenYpix - 1): xylim,...
            -xylim: 2 * xylim / (screenYpix - 1): xylim);
        at = atan2(y, x);
        checks = ((1 + sign(sin(at * LO.tcycles) + eps)...
            .* sign(sin(sqrt(x.^2 + y.^2)))) / 2) * (white - black) + black;
        circle = x.^2 + y.^2 <= xylim^2;
        checks = circle .* checks + grey * ~circle;

        LO.checks = checks;
        LO.radialCheckerboardTexture = Screen('MakeTexture', win, checks);    
    end

	%get the number of conditions
    n_speeds = length(LO.expansion_speeds);
    n_colors = size(LO.colors,1);
    n_alphas = length(LO.alphas);
    
    %number of conditions
    cond_num = n_speeds*n_colors*n_alphas;
    %assemble a matrix with all the parameters to be varied within a trial
    LO.paramorder = zeros(LO.n_reps,cond_num,3);
    
    for r = 1 : LO.n_reps
        switch Param.seqmode
            case 'random'
                [x,y,z] = ind2sub([n_speeds,n_colors,n_alphas],randperm(cond_num));
                LO.paramorder(r,:,:) = permute(cat(1,x,y,z),[3 2 1]);
            case 'sequential'
                [x,y,z] = ind2sub([n_speeds,n_colors,n_alphas],1:cond_num);
                LO.paramorder(r,:,:) = cat(2,x,y,z);
        end
    end
    
    %define the expansion speeds in terms of pix/s
    LO.seqspeeds = pixperdeg.*LO.expansion_speeds(LO.paramorder(:,:,1)).*Param.ifi;
    LO.seqcolors = LO.colors(LO.paramorder(:,:,2)',:);
    LO.seqalphas = LO.alphas(LO.paramorder(:,:,3));
    LO.seqtimes = LO.time_perstimulus(LO.paramorder(:,:,1));
    Param.stimSeq(end+1,1) = { repmat(LO.stim_id, 1,n_colors*n_speeds*n_alphas*LO.n_reps) };
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
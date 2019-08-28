	
% set some DriftingGratings parameters, saved in the structure 'DG'.
function [ DG, Param ] = SetDriftingGratings( win, DG, Param )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

    % Amplitude of the grating in units of absolute display intensity range: A
	% setting of 0.5 means that the grating will extend over a range from -0.5
	% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
	% displayable range. As we select a background color and offset for the
	% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
	% will extend the sinewaves values from 0 = total black in the minima of
	% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
	% than 0.5 don't make sense, as parts of the grating would lie outside the
	% displayable range for your computers displays:
    if ~isfield(DG,'amplitude') || isempty(DG.amplitude)
        DG.amplitude = 0.5;
    end
    
    pixperdeg = Param.pixperdeg(1);
    DG.freq_cppx = DG.spacfreq/pixperdeg;
    if Param.Dichoptic
		if strcmp(inputname(2),'DG2') || strcmp(inputname(2),'DG4')
			pixperdeg = Param.pixperdeg(2);
			DG.freq_cppx = DG.spacfreq/pixperdeg;
        end
    end
	% Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
% 	DG.phase = 0; % I entered it already in the initial parameters
	% Compute increment of phase shift per redraw:
	DG.phaseincrement = (DG.cyclespersecond * 360) * Param.ifi;   % computes increment of degrees per ifi
	% Build a procedural sine grating texture for a grating with a support of
	% gratingsize(1) x gratingsize(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
	% you need to have the mex file CreateProceduralGratingMod and the GSGL
	% shaders.
	DG.gratingRect = [Param.screenRect(1), Param.screenRect(2), Param.screenRect(3), Param.screenRect(4)];
    DG.gratingsize = [Param.screenRect(3), Param.screenRect(4)]; % x,y
    
    
    if DG.drawMask == 1
        DG = SetMask(eye, DG, Param, pixperdeg, win);
% % % 		DG.maskSize_px   = DG.maskSize_deg * pixperdeg ;
% % % 		DG.maskRadius_px = round(DG.maskSize_px/2);
% % %         DG.maskSize_px   = 2*DG.maskRadius_px +1 ;
% % %         
% % %         tmasksize = 2*DG.maskSize_px;
% % %         % calculate offset of grating in pixels:
% % %         DG.maskCenter_offset_px = DG.maskCenter_offset .* [Param.screenRect(3), Param.screenRect(4)] ;
% % %         
% % %         % Important: gratingsize and gratingRect have to be the same size,
% % %         % otherwise the grating will be rescaled
% % %         DG.gratingsize = [tmasksize, tmasksize];
% % %         
% % % %       % Gaussian mask
% % % % 
% % % % 		% Create a single gaussian transparency mask and store it to a texture:
% % % % 		% The mask must have the same size as the visible size of the grating
% % % % 		% to fully cover it. Here we must define it in 2 dimensions and can't
% % % % 		% get easily away with one single row of pixels.
% % % % 		%
% % % % 		% We create a  two-layer texture: One unused luminance channel which we
% % % % 		% just fill with the same color as the background color of the screen
% % % % 		% 'gray'. The transparency (aka alpha) channel is filled with a
% % % % 		% gaussian (exp()) aperture mask:
% % %         % mask = ones(2*DG.maskRadius_px+1, 2*DG.maskRadius_px+1, 2) * DG.BackgroundLuminance;
% % %         % [x,y] = meshgrid(-1*DG.maskRadius_px:1*DG.maskRadius_px, -1*DG.maskRadius_px:1*DG.maskRadius_px);
% % % % 
% % %         % sigma = DG.maskRadius_px * 0.341 ;
% % %         % gaussmask = (1 - exp( -(x.^2 +y.^2)/(2*sigma^2)) *1 );
% % %         % gaussmask(gaussmask<0) = 0;
% % %         % mask(:,:,2) = 255 * gaussmask;
% % %         % [DG.gratingRect,dh,dv] = CenterRect([0 0 DG.maskSize_px DG.maskSize_px], Param.screenRect);
% % %         
% % %         % Smoothed circular aperture mask
% % %         
% % %         % Create transparency mask (enlarged because we will smooth the
% % %         % edge of the actual circular aperture):
% % %         mask = ones(tmasksize, tmasksize, 2) * DG.BackgroundLuminance;
% % %         % Circle with radius DG.maskRadius_px centered at
% % %         % (tmasksize/2,tmasksize/2) in image tmasksizeXtmasksize:
% % %         [rr, cc] = meshgrid(1:tmasksize);
% % %         C = sqrt((rr-tmasksize/2).^2 + (cc-tmasksize/2).^2) <= DG.maskRadius_px ;
% % %         Csm = imfilter(double(C), fspecial('average',[31 31]));
% % %         % Final transparency mask (0=transparent, 255=oapque):
% % %         mask(:,:,2) = 255 * (1 - Csm);
% % % %         figure; imagesc(mask(:,:,2));
% % % 
% % % 		DG.masktex = Screen('MakeTexture', win, mask);
% % % 		
% % % 		% Definition of the drawn rectangle on the screen, centered on the
% % % 		% screen center plus maskCenter_offset:
% % %         [DG.gratingRect,dh,dv] = CenterRect([0 0 tmasksize tmasksize], Param.screenRect);
% % %         DG.gratingRect = DG.gratingRect + [DG.maskCenter_offset_px(1) DG.maskCenter_offset_px(2) DG.maskCenter_offset_px(1) DG.maskCenter_offset_px(2)];
% % % 	
    end

    %create the grating textures
	if DG.sinwave
		DG.gratingtex = CreateProceduralGratingMod(win, 'sin', DG.gratingsize(1), DG.gratingsize(2), [DG.amplitude DG.amplitude DG.amplitude 0], [DG.amplitude DG.amplitude DG.amplitude 0]);
        % default gratingtex, with this you cannot
        % set square but only sine. On the other hand you don't need extra
        % mex files and GSGLShaders:
%         DG.gratingtex = CreateProceduralSineGrating(win, DG.gratingsize(1), DG.gratingsize(2), [0.5 0.5 0.5 0]);
	else
		DG.gratingtex = CreateProceduralGratingMod(win, 'square', DG.gratingsize(1), DG.gratingsize(2), [DG.amplitude DG.amplitude DG.amplitude DG.amplitude], [DG.amplitude DG.amplitude DG.amplitude 0]);
    end
	
    %get the number of conditions
    n_angles = DG.directions;
    n_spatial = length(DG.freq_cppx);
    n_temporal = length(DG.phaseincrement);
    cond_num = n_angles*n_spatial*n_temporal;
    %assemble a matrix with all the parameters to be varied within a trial
    DG.paramorder = zeros(DG.n_reps,cond_num,3);
    
	%Randomize the order of angles (conditions)
    for r = 1 : DG.n_reps
        switch Param.seqmode
            case 'random'
                [x,y,z] = ind2sub([n_angles,n_spatial,n_temporal],randperm(cond_num));
                DG.paramorder(r,:,:) = permute(cat(1,x,y,z),[3 2 1]);
%                     DG.seqdirections(r,:) = randperm(DG.directions);
            case 'sequential'
                [x,y,z] = ind2sub([n_angles,n_spatial,n_temporal],1:cond_num);
                DG.paramorder(r,:,:) = permute(cat(1,x,y,z),[3 2 1]);
%                 DG.seqdirections(r,:) = [1:DG.directions];  % to have sequential sequence of directions, not random.
        end
    end
    %assemble the actual sequence of angles to feed to PTB
%     if ~isempty(DG.angles) && numel(DG.angles)==DG.directions
%         DG.seqangles = repmat(DG.angles+90, DG.n_reps,1) + DG.offset_rot_deg ; % +90 to convert from cartesian to PTB angles
%         angles_cartesian = DG.angles;
%     else
        DG.seqangles = [DG.paramorder(:,:,1)' * 360/DG.directions];
        DG.seqangles(DG.seqangles == 360) = 0;
        DG.seqangles = DG.seqangles + DG.offset_rot_deg ; 
        % list of angles given in cartesian coordinates (0deg=upward,
        % increasing clockwise):
        angles_cartesian = [1:DG.directions]*360/DG.directions-90;
%     end
    DG.angles_cartesian = mod(angles_cartesian, 360);
    DG.seq_spatial = DG.freq_cppx(DG.paramorder(:,:,2)');
    DG.seq_temporal = DG.phaseincrement(DG.paramorder(:,:,3)');
    
    
	Param.stimSeq(end+1,1) = { repmat(DG.stim_id, 1,cond_num*DG.n_reps) };
    DG.cnt=1;
    
if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
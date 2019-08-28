	
% set some DriftingGratingsDisparity parameters, saved in the structure 'DGD'.
function [ DG, Param ] = SetDriftingGratingInDepth( win, DG, Param )
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
    DG.freq_cppx(1) = DG.spacfreq/Param.pixperdeg(1);
    DG.freq_cppx(2) = DG.spacfreq/Param.pixperdeg(2);
%     DG.freq_cppx(2) = DG.freq_cppx(1);

    DG.nSpeeds = numel(DG.cyclespersecond)*2;
    DG.nConditions_perOri = DG.nSpeeds.^2;
    DG.nConditions_tot    = DG.nConditions_perOri * DG.orientations;
    
    DG.speeds = [-fliplr(DG.cyclespersecond) DG.cyclespersecond];
    DG.speedCombinations = AllCombosOf2Vectors(DG.speeds,DG.speeds);
    
	% Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
% 	DG.phase = 0; % I entered it already in the initial parameters
	% Compute increment of phase shift per redraw:
	DG.phaseincrement = (DG.speedCombinations * 360) * Param.ifi;   % computes increment of degrees per ifi
	% Build a procedural sine grating texture for a grating with a support of
	% gratingsize(1) x gratingsize(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
	% you need to have the mex file CreateProceduralGratingMod and the GSGL
	% shaders.
	DG.gratingRect = [Param.screenRect(1), Param.screenRect(2), Param.screenRect(3), Param.screenRect(4)];
    DG.gratingsize = [Param.screenRect(3), Param.screenRect(4)]; % x,y
    
    if DG.drawMask == 1
        eye = 1;
        DG = SetMask(eye, DG, Param, Param.pixperdeg(eye), win);
        eye = 2;
        DG = SetMask(eye, DG, Param, Param.pixperdeg(eye), win);
    end
    
    if DG.sinwave
		DG.gratingtex = CreateProceduralGratingMod(win(1), 'sin', DG.gratingsize(1), DG.gratingsize(2), [DG.amplitude DG.amplitude DG.amplitude 0], [DG.amplitude DG.amplitude DG.amplitude 0]);
			% gratingtex = CreateProceduralSineGrating(win,  res(1), res(2),
			% [0.5 0.5 0.5 0.0]);  % default gratingtex, with this you cannot
			% set square but only sine. On the other hand you don't need extra
			% mex files and GSGLShaders.
	else
		DG.gratingtex = CreateProceduralGratingMod(win(1), 'square', DG.gratingsize(1), DG.gratingsize(2), [DG.amplitude DG.amplitude DG.amplitude 0.0], [DG.amplitude DG.amplitude DG.amplitude 0]);
    end
    
    DG.seqorientations = repmat(1:DG.orientations, DG.nConditions_perOri,1);
    DG.seqconditions = repmat([1:DG.nConditions_perOri]', 1,DG.orientations);
	for r = 1 : DG.n_reps
        switch Param.seqmode
			case 'random'
				DG.stimseq(r,:) = randperm(DG.orientations*DG.nConditions_perOri);
			case 'sequential'
				DG.stimseq(r,:) = 1 : DG.orientations*DG.nConditions_perOri; 
        end
    end
    
    % Phase of left and right screen at t=0:
    DG.InitialPhase = round((359-0)*rand(DG.nConditions_tot, DG.n_reps, 2) + 0);
    
    if ~isempty(DG.angles) && numel(DG.angles)==DG.orientations
        DG.seqangles = repmat(DG.angles+90, DG.nConditions_perOri,1); % +90 to convert from cartesian to PTB angles
        angles_cartesian = DG.angles;
    else
        DG.seqangles = [DG.seqorientations * 180/DG.orientations] ;  % random sequence of orientations (e.g. random sequence of 30:30:360)
        % list of angles given in cartesian coordinates (0deg=upward,
        % increasing clockwise):
        angles_cartesian = [1:DG.orientations]*180/DG.orientations-90;
    end
    DG.seqangles(DG.seqangles == 360) = 0;
    DG.angles_cartesian = mod(angles_cartesian, 180);

	Param.stimSeq(end+1,1) = { repmat(DG.stim_id, 1,DG.n_reps*DG.nConditions_tot) };
    DG.cnt=1;
    
if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end

end

function d = AllCombosOf2Vectors(a,b)

[A,B] = meshgrid(a,b);
c=cat(2,A',B');
d=reshape(c,[],2);

end
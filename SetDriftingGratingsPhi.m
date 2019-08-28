	
% set some DriftingGratings parameters, saved in the structure 'DG'.
function [ DGP, Param ] = SetDriftingGratingsPhi( win, DGP, Param )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

    DGP.patch_time_fr = round(DGP.patch_time/Param.ifi);
    
    % Amplitude of the grating in units of absolute display intensity range: A
	% setting of 0.5 means that the grating will extend over a range from -0.5
	% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
	% displayable range. As we select a background color and offset for the
	% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
	% will extend the sinewaves values from 0 = total black in the minima of
	% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
	% than 0.5 don't make sense, as parts of the grating would lie outside the
	% displayable range for your computers displays:
	DGP.amplitude = 0.5;
    if Param.Dichoptic && strcmp(inputname(2),'DGP2')
        DGP.freq_cppx = DGP.spacfreq/Param.pixperdeg(2);
    else
        DGP.freq_cppx = DGP.spacfreq/Param.pixperdeg(1);
    end
	% Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
% 	DGP.phase = 0; % I entered it already in the initial parameters
% 	% Compute increment of phase shift per redraw:
% 	DGP.phaseincrement = (DGP.cyclespersecond * 360) * Param.ifi;   % computes increment of degrees per ifi
	% Build a procedural sine grating texture for a grating with a support of
	% gratingsize(1) x gratingsize(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
	% you need to have the mex file CreateProceduralGratingMod and the GSGL
	% shaders.
	DGP.gratingsize = [Param.screenRect(1), Param.screenRect(2), Param.screenRect(3), Param.screenRect(4)];
	if DGP.sinwave
		DGP.gratingtex = CreateProceduralGratingMod(win, 'sin', DGP.gratingsize(3), DGP.gratingsize(4), [0.5 0.5 0.5 0], [0.5 0.5 0.5 0]);
			% gratingtex = CreateProceduralSineGrating(win,  res(1), res(2),
			% [0.5 0.5 0.5 0.0]);  % default gratingtex, with this you cannot
			% set square but only sine. On the other hand you don't need extra
			% mex files and GSGLShaders.
	else
		DGP.gratingtex = CreateProceduralGratingMod(win, 'square', DGP.gratingsize(3), DGP.gratingsize(4), [0.5 0.5 0.5 0.0], [0.5 0.5 0.5 0]);
    end
    
    seqdirections = [1:DGP.directions];
    DGP.seqdirections = repmat(seqdirections, length(DGP.phase),1);
    DGP.seqphases = repmat([1:length(DGP.phase)]', 1,DGP.directions);
	for r = 1 : DGP.n_reps
		switch Param.seqmode
			case 'random'
                repeat = true;
                while repeat
                    stimseqtmp = randperm(DGP.directions*length(DGP.phase));
                    dseq = abs(diff(stimseqtmp));
                    % if there are two consecutive stimuli of same
                    % direction or orientation, try again:
                    if any(dseq/DGP.directions/2==fix(dseq/DGP.directions/2))
                        repeat = true;
                    else
                        repeat = false;
                        DGP.stimseq(r,:) = stimseqtmp;
                    end
                end
			case 'sequential'
				DGP.stimseq(r,:) = [1:DGP.directions*length(DGP.phase)];  % to have sequential sequence of directions, not random.
		end
    end

    DGP.seqangles = [DGP.seqdirections * 360/DGP.directions] + DGP.offset_rot_deg ;  % random sequence of directions (e.g. random sequence of 30:30:360)
%     if ~isempty( strfind(inputname(2), 'SG') )
        DGP.seqangles(DGP.seqangles==360) = 0;
%         DGP.seqangles(DGP.seqangles>=180) = DGP.seqangles(DGP.seqangles>=180)-180;
%     end
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    angles_cartesian = seqdirections*360/DGP.directions-90;
    DGP.angles_cartesian = mod(angles_cartesian, 360);
	Param.stimSeq(end+1,1) = { repmat(DGP.stim_id, 1,DGP.directions*length(DGP.phase)*DGP.n_reps) };
    DGP.cnt=1;
    
if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
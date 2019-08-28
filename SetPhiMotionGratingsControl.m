	
% set some DriftingGratings parameters, saved in the structure 'DG'.
function [ PMG, Param ] = SetPhiMotionGratingsControl( win, PMG, Param, screen_id )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
    
    PMG.patch_time    = 1/PMG.cyclespersecond/4; % (0.125s for TF=2Hz)
    PMG.patch_time_fr = round(PMG.patch_time/Param.ifi);
    PMG.frame_interpatch = round(PMG.interpatch_time/Param.ifi);
    % Amplitude of the grating in units of absolute display intensity range: A
	% setting of 0.5 means that the grating will extend over a range from -0.5
	% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
	% displayable range. As we select a background color and offset for the
	% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
	% will extend the sinewaves values from 0 = total black in the minima of
	% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
	% than 0.5 don't make sense, as parts of the grating would lie outside the
	% displayable range for your computers displays:
	PMG.amplitude = 0.5;
    if screen_id==1
        PMG.freq_cppx(1) = PMG.spacfreq/Param.pixperdeg(1);
        % this is to bring the white bar (which is always at the beginning of
        % the texture) to the center of the screen:
        PMG.phase1 = PMG.phase - round( 0.5*Param.screenRect(3)*PMG.freq_cppx(1)*360 ) + 0.25*360;
        % the following two lines are to set the phase shift to bring the white
        % bar to the center of the RF that you set with PMG.offset_phase1_perc:
        offset_phase1_px = PMG.offset_phase1_perc .* [Param.screenRect(3) Param.screenRect(4)];
        PMG.offset_phase1_deg = - round( offset_phase1_px * PMG.freq_cppx(1) * 360 );
    %     PMG.offset_angle1 = atan( offset_phase1_px(2)/offset_phase1_px(1) );
        nr_phases = length(PMG.phase1);
    end
    if screen_id==2
        PMG.freq_cppx(2) = PMG.spacfreq/Param.pixperdeg(2);
        PMG.phase2 = PMG.phase - round( 0.5*Param.screenRect(3)*PMG.freq_cppx(2)*360 ) + 0.25*360 + PMG.phaseincrement ;
        offset_phase2_px = PMG.offset_phase2_perc .* [Param.screenRect(3) Param.screenRect(4)];
        PMG.offset_phase2_deg = - round( offset_phase2_px * PMG.freq_cppx(2) * 360 );
%         PMG.offset_angle2 = atan( offset_phase2_px(2)/offset_phase2_px(1) );
        nr_phases = length(PMG.phase2);
    end
	% Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
% 	PMG.phase = 0; % I entered it already in the initial parameters
% 	% Compute increment of phase shift per redraw:
% 	PMG.phaseincrement = (PMG.cyclespersecond * 360) * Param.ifi;   % computes increment of degrees per ifi
	% Build a procedural sine grating texture for a grating with a support of
	% gratingsize(1) x gratingsize(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
	% you need to have the mex file CreateProceduralGratingMod and the GSGL
	% shaders.
	PMG.gratingsize = [Param.screenRect(1), Param.screenRect(2), Param.screenRect(3), Param.screenRect(4)];
    if PMG.sinwave
		PMG.gratingtex = CreateProceduralGratingMod(win, 'sin', PMG.gratingsize(3), PMG.gratingsize(4), [0.5 0.5 0.5 0], [0.5 0.5 0.5 0]);
			% gratingtex = CreateProceduralSineGrating(win,  res(1), res(2),
			% [0.5 0.5 0.5 0.0]);  % default gratingtex, with this you cannot
			% set square but only sine. On the other hand you don't need extra
			% mex files and GSGLShaders.
	else
		PMG.gratingtex = CreateProceduralGratingMod(win, 'square', PMG.gratingsize(3), PMG.gratingsize(4), [0.5 0.5 0.5 0.0], [0.5 0.5 0.5 0]);
    end    
    
    seqdirections = [1:PMG.directions];
    PMG.seqdirections = repmat(seqdirections, nr_phases,1);
    PMG.seqphases = repmat([1:nr_phases]', 1,PMG.directions);
	for r = 1 : PMG.n_reps
		switch Param.seqmode
			case 'random'
				PMG.stimseq(r,:) = randperm(PMG.directions*nr_phases);
			case 'sequential'
				PMG.stimseq(r,:) = [1:PMG.directions*nr_phases];  % to have sequential sequence of directions, not random.
		end
	end

    PMG.seqangles = [PMG.seqdirections * 360/PMG.directions];  % random sequence of directions (e.g. random sequence of 30:30:360)
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    angles_cartesian = [1:PMG.directions]*360/PMG.directions-90;
    PMG.angles_cartesian = mod(angles_cartesian, 360);
    
	Param.stimSeq(end+1,1) = { repmat(PMG.stim_id, 1,PMG.directions*nr_phases*PMG.n_reps) };
    PMG.cnt=1;
    
if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
	
% set some DriftingGratings parameters, saved in the structure 'DG'.
function [ PMG, Param ] = SetPhiMotionGratings( win, PMG, Param )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
    
%     PMG.patch_time    = 1/PMG.cyclespersecond/4; % (0.125s for TF=2Hz)
%     PMG.patch_time_fr = round(PMG.patch_time/Param.ifi);
%     PMG.frame_interpatch = round(PMG.interpatch_time/Param.ifi);
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
    PMG.freq_cppx(1) = PMG.spacfreq/Param.pixperdeg(1);
    PMG.freq_cppx(2) = PMG.spacfreq/Param.pixperdeg(2);
% %     % this is to bring the white bar (which is always at the beginning of
% %     % the texture) to the center of the screen:
% %     PMG.phase1 = PMG.phase - round( 0.5*Param.screenRect(3)*PMG.freq_cppx(1)*360 ) + 0.25*360;
%     PMG.phase1 = PMG.phase;
    offset_toscreencenter1_deg = - round( 0.5*Param.screenRect(3)*PMG.freq_cppx(1)*360 ) + 0.25*360;
    offset_toscreencenter2_deg = - round( 0.5*Param.screenRect(3)*PMG.freq_cppx(2)*360 ) + 0.25*360;
    % The following two lines are to set the phase shift to bring the white
    % bar to the center of the RF with phase=0 (you set RF center with PMG.offset_phase1_perc):
    offset_phase1_px = PMG.offset_phase1_perc .* [Param.screenRect(3) Param.screenRect(4)];
    offset_phase2_px = PMG.offset_phase2_perc .* [Param.screenRect(3) Param.screenRect(4)];
    PMG.offset_phase1_deg = - round( offset_phase1_px * PMG.freq_cppx(1) * 360 ) + offset_toscreencenter1_deg;
    PMG.offset_phase2_deg = - round( offset_phase2_px * PMG.freq_cppx(2) * 360 ) + offset_toscreencenter2_deg;
%     PMG.offset_angle1 = atan( offset_phase1_px(2)/offset_phase1_px(1) );
%     if Param.Dichoptic

% %         PMG.phase2 = PMG.phase - round( 0.5*Param.screenRect(3)*PMG.freq_cppx(2)*360 ) + 0.25*360 + PMG.phaseincrement ;
%         PMG.phase2 = PMG.phase;
        
        
%         PMG.offset_angle2 = atan( offset_phase2_px(2)/offset_phase2_px(1) );
%     end
    
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
%         PMG.gratingtex = CreateProceduralGratingMod(win, 'sin', PMG.gratingsize(3), PMG.gratingsize(4), [0.5 0.5 0.5 0.0], [0.5 0.5 0.5 0]);
		PMG.gratingtex = CreateProceduralGratingMod2(win, 'sin', PMG.gratingsize(3), PMG.gratingsize(4), [0.5 0.5 0.5 0]);
			% gratingtex = CreateProceduralSineGrating(win,  res(1), res(2),
			% [0.5 0.5 0.5 0.0]);  % default gratingtex, with this you cannot
			% set square but only sine. On the other hand you don't need extra
			% mex files and GSGLShaders.
    else
%         PMG.gratingtex = CreateProceduralGratingMod(win, 'square', PMG.gratingsize(3), PMG.gratingsize(4), [0.5 0.5 0.5 0.0], [0.5 0.5 0.5 0]);
		PMG.gratingtex = CreateProceduralGratingMod2(win, 'square', PMG.gratingsize(3), PMG.gratingsize(4), [0.5 0.5 0.5 0.0]);
    end    
    
    nCombos = length(PMG.phases)/2*length(PMG.phases); % nr. combos per orientations (includes both directions per ori.):
    PMG.nCombos = nCombos;
    [A,B]=meshgrid( PMG.phases(1:end/2) , PMG.phases );
    seqphases = reshape( cat(2,A,B) ,[],2 );
    seqphases = repmat(seqphases,PMG.orientations,1);
    PMG.seqphases = seqphases; % size (nCombos,2)
    
    seqorientations = [1:PMG.orientations];
    seqorientations = repmat(seqorientations, nCombos,1);
    PMG.seqorientations = seqorientations;
%     PMG.seqphases = repmat([1:length(PMG.phase1)]', 1,PMG.directions);
	for r = 1 : PMG.n_reps %#ok<*ALIGN>
		switch Param.seqmode
			case 'random'
				stimseq(r,:) = randperm(PMG.orientations*nCombos);  %#ok<*AGROW>
			case 'sequential'
				stimseq(r,:) = [1:PMG.orientations*nCombos];  % to have sequential sequence of directions, not random.
        end
    end
    PMG.stimseq = stimseq'; % (check that stimseq columns => reps)
    
    if ~isempty(PMG.angles) && numel(PMG.angles)==PMG.orientations
        seqangles = repmat(PMG.angles+90, nCombos,1) ; % +90 to convert from cartesian to PTB angles
        angles_cartesian = PMG.angles;
    else
        seqangles = [PMG.seqorientations * 180/PMG.orientations];  % random sequence of directions (e.g. random sequence of 30:30:360)
        % list of angles given in cartesian coordinates (0deg=upward,
        % increasing clockwise):
        angles_cartesian = [1:PMG.orientations]*180/PMG.orientations-90;
    end
    seqangles(seqangles==180) = 0;
    PMG.seqangles = seqangles; % (check that seqangles columns => different orientations)
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    PMG.angles_cartesian = mod(angles_cartesian, 180);

	Param.stimSeq(end+1,1) = { repmat(PMG.stim_id, 1,PMG.orientations*nCombos*PMG.n_reps) };
    PMG.cnt=1;
    
if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
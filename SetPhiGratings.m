	
% set some PhiGratings parameters, saved in the structure 'PG'.
function [ PG, Param ] = SetPhiGratings( win, PG, Param )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
    

%     PG.patch_time    = 1/PG.cyclespersecond/4; % (0.125s for TF=2Hz)
    PG.patch_time_fr = round(PG.patch_time/Param.ifi);
    PG.frame_interpatch = round(PG.interpatch_time/Param.ifi);
    % Amplitude of the grating in units of absolute display intensity range: A
	% setting of 0.5 means that the grating will extend over a range from -0.5
	% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
	% displayable range. As we select a background color and offset for the
	% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
	% will extend the sinewaves values from 0 = total black in the minima of
	% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
	% than 0.5 don't make sense, as parts of the grating would lie outside the
	% displayable range for your computers displays:
	PG.amplitude = 0.5;
    PG.freq_cppx(1) = PG.spacfreq/Param.pixperdeg(1);
    PG.freq_cppx(2) = PG.spacfreq/Param.pixperdeg(2);

	% Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
% 	PG.phase = 0; % I entered it already in the initial parameters
% 	% Compute increment of phase shift per redraw:
% 	PG.phaseincrement = (PG.cyclespersecond * 360) * Param.ifi;   % computes increment of degrees per ifi
	% Build a procedural sine grating texture for a grating with a support of
	% gratingsize(1) x gratingsize(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
	% you need to have the mex file CreateProceduralGratingMod and the GSGL
	% shaders.
	PG.gratingsize = [Param.screenRect(1), Param.screenRect(2), Param.screenRect(3), Param.screenRect(4)];
    if PG.sinwave
		PG.gratingtex = CreateProceduralGratingMod(win, 'sin', PG.gratingsize(3), PG.gratingsize(4), [0.5 0.5 0.5 0], [0.5 0.5 0.5 0]);
			% gratingtex = CreateProceduralSineGrating(win,  res(1), res(2),
			% [0.5 0.5 0.5 0.0]);  % default gratingtex, with this you cannot
			% set square but only sine. On the other hand you don't need extra
			% mex files and GSGLShaders.
	else
		PG.gratingtex = CreateProceduralGratingMod(win, 'square', PG.gratingsize(3), PG.gratingsize(4), [0.5 0.5 0.5 0.0], [0.5 0.5 0.5 0]);
    end    
    
    
    
    % nr of combinations per orientation (includes both directions):
    PG.nCombinations = factorial(PG.nbr_phases)/factorial(2); % =12 when nbr_phases=4.
    
% PG.seqphases contains the sequence of phases 1 1; 2 2; 3 3; 4 4; 1 2; 1 3...
% until 4 3. This sequence is first run dichoptically for screensequence 1 2, then the
% same for screensequence 2 1. This block is then repeated for each
% orientation. This last block is then run monocularly, with screensequence
% 1 1 and then 2 2.
% 
PG.seqphases=[]; seqphasesdichtmp=[]; seqphasesdich=[]; rowstokeep=[];
[A,B]=meshgrid( 1:length(PG.phases) , 1:length(PG.phases) );
seqphasestmp = reshape( cat(2,A,B) ,[],2 ); % sequences of dichoptic phi motion bars and "pseudostatic" bars
for ii=1:size(seqphasestmp,1)
    if seqphasestmp(ii,1)~=seqphasestmp(ii,2); rowstokeep = [rowstokeep,ii]; end;
end; seqphasestmp = seqphasestmp(rowstokeep,:);
% if Param.Dichoptic
% %     additionalStaticBars = repmat([1,1; 2,2; 3,3] ,3,1);
    PseudoStaticBars = 1:length(PG.phases); PseudoStaticBars = [PseudoStaticBars;PseudoStaticBars]';
%     seqphasestmp=cat(1,seqphasestmp,additionalStaticBars);
% end
seqphasesdichtmp = [ PseudoStaticBars; seqphasestmp];
seqphasesdichtmp = [seqphasesdichtmp ; seqphasesdichtmp]; % duplicate for screensequence [1 2] and [2 1]

seqphasesdich = repmat(seqphasesdichtmp, PG.directions/2, 1) ;

seqphasesmono = seqphasesdich;
PG.seqphases = [seqphasesdich; seqphasesmono];

    screensequencedich = [ repmat([1 2], length(seqphasesdichtmp)/2,1); repmat([2 1], length(seqphasesdichtmp)/2,1)];
    screensequencemono = [screensequencedich(:,1),screensequencedich(:,1)];
    % the following line is to have only a single flash for monoc stimuli
    % with same phase (=static stimuli):
    %   screensequencemono([1:length(PG.phases), length(PG.phases)^2+1:length(PG.phases)^2+length(PG.phases)],2)=0;
    screensequencedich = repmat( screensequencedich, PG.directions/2,1);
    screensequencemono = repmat( screensequencemono, PG.directions/2,1);

    PG.screensequence = [ screensequencedich; screensequencemono];

    seqorientations = [];
    for d = 1 : PG.directions/2
        seqorientations = cat( 1, seqorientations, repmat(d, (PG.nbr_phases^2)*2,1) );
    end
    % replicate twice for dichoptic and monoc:
    PG.seqorientations = [seqorientations; seqorientations];
    
for r = 1:PG.n_reps
    switch Param.seqmode
        case 'sequential'
            PG.stimseq(r,:) = 1 : PG.nbr_stimuli ;
        case 'random'
            % this 'random' prevents that the same orientation (regardless
            % of the specific pattern) is presented twice (or more times)
            % in a row:
            stimseq = [];
            for d = 1 : PG.directions/2 *2 % /2 because we want the orientations, *2 because we have dichoptic and then monoc
                randstimseq_perOri(d,:) = (d-1)*(PG.nbr_phases^2)*2 + randperm( (PG.nbr_phases^2)*2 ); %#ok<AGROW>
            end
            for i = 1:size(randstimseq_perOri,2)
                repeat = true;
                while repeat
                    randtmp(i,:) = randperm( PG.directions/2 *2 ); %#ok<AGROW>
                    randtmp_diff = abs(diff(randtmp(i,:)));
                    if any(randtmp_diff == (PG.directions/2))
                        repeat=true;
                    else
                        if i>1 && ...
                            (randtmp(i,1) == randtmp(i-1,end) || ...
                                abs(randtmp(i,1) - randtmp(i-1,end))==(PG.directions/2) )
                            repeat=true;
                        else
                            repeat= false;
                        end
                    end
                        
                end
                stimseq_column = randstimseq_perOri(randtmp(i,:),i);
                stimseq = cat(2, stimseq, stimseq_column);
            end
            PG.stimseq(r,:) = reshape( stimseq, [],1 ) ;
    end     
end

    PG.seqangles = [PG.seqorientations * 360/PG.directions];  %sequence of directions (e.g. random sequence of 30:30:360)
    PG.seqangles(PG.seqangles==180) = 0 ;
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    angles_cartesian = [1:(PG.directions/2)]*360/PG.directions-90;
    PG.angles_cartesian = mod(angles_cartesian, 360);
    
Param.stimSeq(end+1,1) = { repmat(PG.stim_id, 1,PG.nbr_stimuli*PG.n_reps) };
PG.cnt=1; 

if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
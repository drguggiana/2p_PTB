
function [ RP, Param ] = SetRetinotopicPosition( win, RP, Param )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if sum(RP.margin.X) >= RP.tot_positions(1) ||...
    sum(RP.margin.Y) >= RP.tot_positions(2)
    error('  /!\  Margins are larger than tot_positions: did you set RP.margin and RP.tot_positions correctly??');
end
% set some RetinotopicPosition parameters, saved in the structure 'RP'.

    RP.stimulus_time_fr = round(RP.stimulus_time/Param.ifi);
	RP.amplitude = 0.5;
    if Param.Dichoptic && strcmp(inputname(2),'RP2')
        RP.freq_cppx = RP.spacfreq/Param.pixperdeg(2) ;
    else
        RP.freq_cppx = RP.spacfreq/Param.pixperdeg(1) ;
    end
	RP.phase = 0;
	RP.phaseincrement = (RP.cyclespersecond * 360) * Param.ifi;
	RP.positions = [RP.tot_positions(1)-sum(RP.margin.X) , RP.tot_positions(2)-sum(RP.margin.Y)];
	RP.gratingsize = round([ Param.screenRect(3)/RP.tot_positions(1) , Param.screenRect(4)/RP.tot_positions(2) ]) ;
    for py = 1 : RP.tot_positions(2)
	for px = 1 : RP.tot_positions(1)
		RP.coord_positions(px,py,:) = [ Param.screenRect(1) + (px-1)*RP.gratingsize(1) ...
										Param.screenRect(2) + (py-1)*RP.gratingsize(2) ...
										Param.screenRect(1) + (px)  *RP.gratingsize(1) ...
										Param.screenRect(2) + (py)  *RP.gratingsize(2) ] ;
	end
    end
	RP.seqpositions=[];
    for r = 1:RP.n_reps
		% seqpositionstmp is [1 1; 1 2; 1 3; 1 4; 2 1; ...; 6 4]:
        [A,B]=meshgrid( 1+RP.margin.X(1):RP.tot_positions(1)-RP.margin.X(2) , 1+RP.margin.Y(1):RP.tot_positions(2)-RP.margin.Y(2) );
        seqpositionstmp = reshape( cat(2,A',B') ,[],2 );
		switch Param.seqmode
			case 'random'
				seqpositionstmp=seqpositionstmp(randperm(size(seqpositionstmp,1)),:);
		end
		RP.seqpositions=[ RP.seqpositions ; seqpositionstmp ];
    end
    if RP.sinwave
		RP.gratingtex = CreateProceduralGratingMod(win, 'sin', RP.gratingsize(1), RP.gratingsize(2), [0.5 0.5 0.5 0], [0.5 0.5 0.5 0]);
	else
		RP.gratingtex = CreateProceduralGratingMod(win, 'square', RP.gratingsize(1), RP.gratingsize(2), [0.5 0.5 0.5 0.0], [0.5 0.5 0.5 0]);
    end
    
%     stimSeqtmp3 = [ repmat(RP.stim_id, 1,prod(RP.positions)*RP.n_reps) ];
    Param.stimSeq(end+1,1) = { repmat(RP.stim_id, 1,prod(RP.positions)*RP.n_reps) };
    RP.tempo=cell(1);%prod(RP.positions)*RP.n_reps,1);
    RP.cnt=1;
    
    if exist('scriptName','var')
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
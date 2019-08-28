
function [ RPS, Param ] = SetRetinotopicPositionStripes( win, RPS, Param )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if RPS.tot_positions(1) == 0
    RPS.margin.X = [0 0];
end
if RPS.tot_positions(2) == 0
    RPS.margin.Y = [0 0];
end

if RPS.tot_positions(1) > 0 && RPS.tot_positions(2) > 0
    if sum(RPS.margin.X) >= RPS.tot_positions(1) ||...
        sum(RPS.margin.Y) >= RPS.tot_positions(2)
        error('  /!\  Margins are larger than tot_positions: did you set RPS.margin and RPS.tot_positions correctly??');
    end
end
% set some RetinotopicPosition parameters, saved in the structure 'RPS'.

    RPS.stimulus_time_fr = round(RPS.stimulus_time/Param.ifi);
    if ~isfield(RPS,'amplitude') || isempty(RPS.amplitude)
        RPS.amplitude = 0.5;
    end
    if Param.Dichoptic && strcmp(inputname(2),'RPS2')
        RPS.freq_cppx = RPS.spacfreq/Param.pixperdeg(2) ;
    else
        RPS.freq_cppx = RPS.spacfreq/Param.pixperdeg(1) ;
    end
	RPS.phase = 0;
	RPS.phaseincrement = (RPS.cyclespersecond * 360) * Param.ifi;
	RPS.positions = [RPS.tot_positions(1)-sum(RPS.margin.X) , RPS.tot_positions(2)-sum(RPS.margin.Y)];
    if RPS.OverlapStripes
        RPS.positions = RPS.positions + (RPS.positions-1);
    end
    RPS.positions(RPS.positions<0) = 0;
% 	RPS.gratingsize = round([ Param.screenRect(3)/RPS.tot_positions(1) , Param.screenRect(4)/RPS.tot_positions(2) ]) ;
    RPS.gratingsize = [ Param.screenRect(3) , Param.screenRect(4) ] ;
    pp=0;
% 	for px = 1+RPS.margin.X(1) : RPS.tot_positions(1)-RPS.margin.X(2)
% 		pp=pp+1;
% 		RPS.coord_positions{pp} = [ Param.screenRect(1) + (px-1)*RPS.gratingsize(1) ...
% 									Param.screenRect(2)                             ...
% 									Param.screenRect(1) + (px)  *RPS.gratingsize(1) ...
% 									Param.screenRect(4) ] ;
%     
%         RPS.coord_positions{pp}(RPS.coord_positions{pp}==0) = 1 ;
% 	end
% 	for py = 1+RPS.margin.Y(1) : RPS.tot_positions(2)-RPS.margin.Y(2)
% 		pp=pp+1;
% 		RPS.coord_positions{pp} = [ Param.screenRect(1)                             ...
% 									Param.screenRect(2) + (py-1)*RPS.gratingsize(2) ...
% 									Param.screenRect(3)                             ...
% 									Param.screenRect(2) + (py)  *RPS.gratingsize(2) ] ;
%     end
    
% coordinates for vertical stripes
indx = round( linspace(1,Param.screenRect(3), RPS.tot_positions(1)+1) );
indx_nomarg = indx(1+RPS.margin.X(1) : RPS.tot_positions(1)-RPS.margin.X(2));
widthx_px = mode( diff(indx) );
if RPS.OverlapStripes
    indx2 = round( indx_nomarg(1:end-1) + (indx_nomarg(2:end)-indx_nomarg(1:end-1)) / 2 );
    indx_nomarg_l = sort([indx_nomarg, indx2]);
else
    indx_nomarg_l = indx_nomarg;
end
indx_nomarg_r = indx_nomarg_l + widthx_px;
% coordinates for horizontal stripes
indy = round( linspace(1,Param.screenRect(4), RPS.tot_positions(2)+1) );
indy_nomarg = indy(1+RPS.margin.Y(1) : RPS.tot_positions(2)-RPS.margin.Y(2));
widthy_px = mode( diff(indy) );
if RPS.OverlapStripes
    indy2 = round( indy_nomarg(1:end-1) + (indy_nomarg(2:end)-indy_nomarg(1:end-1)) / 2 );
    indy_nomarg_t = sort([indy_nomarg, indy2]);
else
    indy_nomarg_t = indy_nomarg;
end
indy_nomarg_b = indy_nomarg_t + widthy_px;
%
    mask = ones(RPS.gratingsize(2),RPS.gratingsize(1),2) * RPS.BackgroundLuminance ;
    mask(:,:,2) = 255;
    if RPS.tot_positions(1) > 0 && RPS.positions(1) > 0
        for px = 1 : RPS.positions(1)
            pp = pp+1;
            tmpmask = mask;
            tmpmask(:,indx_nomarg_l(px):indx_nomarg_r(px),2) = 0;
            RPS.maskgrating{pp} = tmpmask;
            RPS.masktex{pp} = Screen('MakeTexture', win, tmpmask);
        end
    end
    if RPS.tot_positions(2) > 0 && RPS.positions(2) > 0
        for py = 1 : RPS.positions(2)
            pp=pp+1;
            tmpmask = mask;
            tmpmask(indy_nomarg_t(py):indy_nomarg_b(py),:,2) = 0;
            RPS.maskgrating{pp} = tmpmask;
            RPS.masktex{pp} = Screen('MakeTexture', win, tmpmask);
        end
    end
    
	RPS.seqpositions=[];
    for r = 1:RPS.n_reps
		% seqpositionstmp is [1; 2; 3;... sum(positions)]:
		switch Param.seqmode
			case 'random'
				seqpositionstmp = randperm(sum(RPS.positions));
			case 'sequential'
				seqpositionstmp = 1 : sum(RPS.positions);
		end
		RPS.seqpositions=[ RPS.seqpositions ; seqpositionstmp ];
    end
    RPS.seqpositions = RPS.seqpositions';
    
    if RPS.sinwave
		RPS.gratingtex = CreateProceduralGratingMod(win, 'sin', RPS.gratingsize(1), RPS.gratingsize(2), [RPS.amplitude RPS.amplitude RPS.amplitude 0], [RPS.amplitude RPS.amplitude RPS.amplitude 0]);
	else
		RPS.gratingtex = CreateProceduralGratingMod(win, 'square', RPS.gratingsize(1), RPS.gratingsize(2), [RPS.amplitude RPS.amplitude RPS.amplitude 0], [RPS.amplitude RPS.amplitude RPS.amplitude 0]);
    end
    
    Param.stimSeq(end+1,1) = { repmat(RPS.stim_id, 1,sum(RPS.positions)*RPS.n_reps) };
    RPS.tempo=cell(1);%sum(RPS.positions)*RPS.n_reps,1);
    RPS.cnt=1;
    
    if exist('scriptName','var')
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
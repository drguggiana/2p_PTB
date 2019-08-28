function [ RBS, Param ] = SetRandomBarStereogram( RBS, Param )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if isempty(RBS.apertureSize)
    RBS.apertureSize = fix([Param.screenRect(3),Param.screenRect(4)]/Param.pixperdeg);
end

RBS.frames_stim     = round(RBS.stimulus_time/Param.ifi); 
RBS.frames_poststim = round(RBS.poststim_time/Param.ifi);
RBS.pattern_time_fr = round(RBS.pattern_time/Param.ifi);

RBS.center = [0,0]; % center of the field of dots (x,y)

disparity_bins = linspace(RBS.disparity_range(1),RBS.disparity_range(2), RBS.n_disparity_intervals);
bin_width = diff(RBS.disparity_range)/(RBS.n_disparity_intervals-1);
bin_ranges_deg = [(disparity_bins-bin_width/2)', (disparity_bins+bin_width/2)'];
% bin_ranges_px = round(bin_ranges * Param.pixperdeg);
RBS.bin_width     = bin_width;
RBS.bin_ranges_deg  = bin_ranges_deg;
% RBS.bin_ranges_px = bin_ranges_px;


RBS.barWidth_px  = round( RBS.barWidth .* Param.pixperdeg(1) );
RBS.barHeight_px = 10000;

% BM.color_list = uint8(BM.color{1});
if     length(RBS.color)==2
    %half nDots with color{1} and half with color{2}:
    colorlist1 = repmat(RBS.color{1}', 1,ceil(RBS.nBars/2));
    colorlist2 = repmat(RBS.color{2}', 1, fix(RBS.nBars/2));
    RBS.color_list = uint8([colorlist1, colorlist2]);
elseif length(RBS.color)==1
    RBS.color_list = uint8(repmat(RBS.color{1}', 1,RBS.nBars));
end


RBS.lifetime_fr = round(RBS.lifetime/Param.ifi);


RBS.seqorientations = repmat(1:RBS.orientations, RBS.n_disparity_intervals,1);
RBS.seqIntervals = repmat([1:RBS.n_disparity_intervals]', 1,RBS.orientations);
for r = 1 : RBS.n_reps
    switch Param.seqmode
        case 'random'
            RBS.stimseq(r,:) = randperm(RBS.orientations*RBS.n_disparity_intervals);
        case 'sequential'
            RBS.stimseq(r,:) = 1 : RBS.orientations*RBS.n_disparity_intervals; 
    end
end

if ~isempty(RBS.angles) && numel(RBS.angles)==RBS.orientations
    RBS.seqangles = repmat(RBS.angles+0, RBS.n_disparity_intervals,1); % +90 to convert from cartesian to PTB angles
    angles_cartesian = RBS.angles;
else
    RBS.seqangles = [RBS.seqorientations * 180/RBS.orientations] ;  % random sequence of orientations (e.g. random sequence of 30:30:360)
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    angles_cartesian = [1:RBS.orientations]*180/RBS.orientations-0;
end
RBS.seqangles(RBS.seqangles == 180) = 0;
RBS.angles_cartesian = mod(angles_cartesian, 180);


Param.stimSeq(end+1,1) = { repmat(RBS.stim_id, 1,RBS.orientations*RBS.n_disparity_intervals*RBS.n_reps) };
RBS.cnt=1;


if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
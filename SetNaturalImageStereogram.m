function [ RDS, Param, img ] = SetNaturalImageStereogram( RDS, Param ) %#ok<STOUT>

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

% load variable img (cell array of nr. 1000 uint8 images):
load([pwd filesep 'img1.mat']);

RDS.frames_stim     = round(RDS.stimulus_time/Param.ifi); 
RDS.frames_poststim = round(RDS.poststim_time/Param.ifi);
RDS.pattern_time_fr = round(RDS.pattern_time/Param.ifi);
RDS.nPatterns = ceil(RDS.stimulus_time / RDS.pattern_time);

disparity_bins = linspace(RDS.disparity_range(1),RDS.disparity_range(2), RDS.n_disparity_intervals);
bin_width = diff(RDS.disparity_range)/(RDS.n_disparity_intervals-1);
bin_ranges_deg = [(disparity_bins-bin_width/2)', (disparity_bins+bin_width/2)'];
% bin_ranges_px = round(bin_ranges * Param.pixperdeg);
RDS.bin_width     = bin_width;
RDS.bin_ranges_deg  = bin_ranges_deg;
% RDS.bin_ranges_px = bin_ranges_px;


RDS.seqorientations = repmat(1:RDS.orientations, RDS.n_disparity_intervals,1);
RDS.seqIntervals = repmat([1:RDS.n_disparity_intervals]', 1,RDS.orientations);
for r = 1 : RDS.n_reps
    switch Param.seqmode
        case 'random'
            RDS.stimseq(r,:) = randperm(RDS.orientations*RDS.n_disparity_intervals);
        case 'sequential'
            RDS.stimseq(r,:) = 1 : RDS.orientations*RDS.n_disparity_intervals; 
    end
end


if ~isempty(RDS.angles) && numel(RDS.angles)==RDS.orientations
    RDS.seqangles = repmat(RDS.angles+0, RDS.n_disparity_intervals,1); % +90 to convert from cartesian to PTB angles
    angles_cartesian = RDS.angles;
else
    RDS.seqangles = [RDS.seqorientations * 180/RDS.orientations] ;  % random sequence of orientations (e.g. random sequence of 30:30:360)
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    angles_cartesian = [1:RDS.orientations]*180/RDS.orientations-0;
end
RDS.seqangles(RDS.seqangles == 180) = 0;
RDS.angles_cartesian = mod(angles_cartesian, 180);


Param.stimSeq(end+1,1) = { repmat(RDS.stim_id, 1,RDS.orientations*RDS.n_disparity_intervals*RDS.n_reps) };
RDS.cnt=1;


if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
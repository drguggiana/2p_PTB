function [ RDS, Param ] = SetRandomDotStereogram( RDS, Param )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if isempty(RDS.apertureSize)
    RDS.apertureSize = fix([Param.screenRect(3),Param.screenRect(4)]/mean(Param.pixperdeg));
end

RDS.frames_stim     = round(RDS.stimulus_time/Param.ifi); 
RDS.frames_poststim = round(RDS.poststim_time/Param.ifi);
RDS.pattern_time_fr = round(RDS.pattern_time/Param.ifi);

disparity_bins = linspace(RDS.disparity_range(1),RDS.disparity_range(2), RDS.n_disparity_intervals);
bin_width = diff(RDS.disparity_range)/(RDS.n_disparity_intervals-1);
bin_ranges_deg = [(disparity_bins-bin_width/2)', (disparity_bins+bin_width/2)'];
% bin_ranges_px = round(bin_ranges * Param.pixperdeg);
RDS.bin_width     = bin_width;
RDS.bin_ranges_deg  = bin_ranges_deg;
% RDS.bin_ranges_px = bin_ranges_px;

RDS.center = [0,0]; % center of the field of dots (x,y)

RDS.size_px = round( RDS.size * Param.pixperdeg );
areaScreen = prod(Param.screenRes);
areaDot = (mean(RDS.size_px)/2)^2*pi;
if isempty(RDS.nDots)
    RDS.nDots = round( areaScreen*RDS.density/areaDot );
end

if     length(RDS.color)==2 && ~isempty(RDS.color{2})
    %half nDots with color{1} and half with color{2}:
    colorlist1 = repmat(RDS.color{1}', 1,ceil(RDS.nDots/2));
    colorlist2 = repmat(RDS.color{2}', 1, fix(RDS.nDots/2));
    RDS.color_list{1} = uint8([colorlist1, colorlist2]);
    if RDS.Anticorrelated == 1
    RDS.color_list{2} = uint8([colorlist2, colorlist1]);
    else
    RDS.color_list{2} = uint8([colorlist1, colorlist2]); 
    end
elseif length(RDS.color)==1 || isempty(RDS.color{2})
    RDS.color_list{1} = uint8(repmat(RDS.color{1}', 1,RDS.nDots));
    RDS.color_list{2} = RDS.color_list{1};
end

% We need to deal with DM moving beyond the edge of the aperture. This
% requrires a couple more lines of code.
% First we'll calculate the left,
% right top and bottom of the aperture (in degrees)
aperture.l = RDS.center(1)-RDS.apertureSize(1)/2;
aperture.r = RDS.center(1)+RDS.apertureSize(1)/2;
aperture.b = RDS.center(2)-RDS.apertureSize(2)/2;
aperture.t = RDS.center(2)+RDS.apertureSize(2)/2;
RDS.aperture = aperture;

RDS.lifetime_fr = round(RDS.lifetime/Param.ifi);
% Each dot will have a integer value 'life' which is how many frames the %
% dot has been going. The starting 'life' of each dot will be a random %
% number between 0 and DM.lifetime-1 so that they don't all 'die' on the
% same frame:
% DM.life = ceil(rand(1,DM.nDots)*DM.lifetime_fr);

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
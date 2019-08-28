function [ BM, Param ] = SetBarMotion( BM, Param )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if isempty(BM.apertureSize)
    BM.apertureSize = fix([Param.screenRect(3),Param.screenRect(4)]/Param.pixperdeg);
end

BM.frames_stim     = round(BM.stimulus_time/Param.ifi); 
BM.frames_poststim = round(BM.poststim_time/Param.ifi);
if isfield(BM,'phi_interval') % for DotMotionPhi
    BM.phi_interval_fr = round(BM.phi_interval/Param.ifi);
end

BM.center = [0,0]; % center of the field of dots (x,y)

BM.barWidth_px  = round( BM.barWidth .* Param.pixperdeg );
BM.barHeight_px = 10000;

% BM.color_list = uint8(BM.color{1});
if     length(BM.color)==2
    %half nDots with color{1} and half with color{2}:
    colorlist1 = repmat(BM.color{1}', 1,ceil(BM.nBars/2));
    colorlist2 = repmat(BM.color{2}', 1, fix(BM.nBars/2));
    BM.color_list = uint8([colorlist1, colorlist2]);
elseif length(BM.color)==1
    BM.color_list = uint8(repmat(BM.color{1}', 1,BM.nBars));
end


BM.lifetime_fr = round(BM.lifetime/Param.ifi);


for r = 1 : BM.n_reps
    switch Param.seqmode
        case 'random'
            BM.seqdirections(r,:) = randperm(BM.directions);
        case 'sequential'
            BM.seqdirections(r,:) = [1:BM.directions];  % to have sequential sequence of directions, not random.
    end
end

BM.seqangles = [BM.seqdirections * 360/BM.directions] + BM.offset_rot_deg ;  % random sequence of directions (e.g. random sequence of 30:30:360)
% list of angles given in cartesian coordinates (0deg=upward,
% increasing clockwise):
angles_cartesian = [1:BM.directions]*360/BM.directions +90;
BM.angles_cartesian = mod(angles_cartesian, 360);

Param.stimSeq(end+1,1) = { repmat(BM.stim_id, 1,BM.directions*BM.n_reps) };
BM.cnt=1;


if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
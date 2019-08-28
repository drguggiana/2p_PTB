function [ DM, Param ] = SetDotMotion( DM, Param )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

if isempty(DM.apertureSize)
    DM.apertureSize = fix([Param.screenRect(3),Param.screenRect(4)]/Param.pixperdeg(1));
end

DM.frames_stim     = round(DM.stimulus_time/Param.ifi); 
DM.frames_poststim = round(DM.poststim_time/Param.ifi);
if isfield(DM,'phi_interval') % for DotMotionPhi
    DM.phi_interval_fr = round(DM.phi_interval/Param.ifi);
end

DM.center = [0,0]; % center of the field of dots (x,y)

DM.size_px = round( DM.size * Param.pixperdeg(1) );
areaScreen = prod(Param.screenRes);
areaDot = (DM.size_px/2)^2*pi;
if isempty(DM.nDots)
    DM.nDots = round( areaScreen*DM.density/areaDot );
end

if     length(DM.color)==2
    %half nDots with color{1} and half with color{2}:
    colorlist1 = repmat(DM.color{1}', 1,ceil(DM.nDots/2));
    colorlist2 = repmat(DM.color{2}', 1, fix(DM.nDots/2));
    DM.color_list = uint8([colorlist1, colorlist2]);
elseif length(DM.color)==1
    DM.color_list = uint8(repmat(DM.color{1}', 1,DM.nDots));
end

% We need to deal with DM moving beyond the edge of the aperture. This
% requrires a couple more lines of code.
% First we'll calculate the left,
% right top and bottom of the aperture (in degrees)
aperture.l = DM.center(1)-DM.apertureSize(1)/2;
aperture.r = DM.center(1)+DM.apertureSize(1)/2;
aperture.b = DM.center(2)-DM.apertureSize(2)/2;
aperture.t = DM.center(2)+DM.apertureSize(2)/2;
DM.aperture = aperture;

DM.lifetime_fr = round(DM.lifetime/Param.ifi);
% Each dot will have a integer value 'life' which is how many frames the %
% dot has been going. The starting 'life' of each dot will be a random %
% number between 0 and DM.lifetime-1 so that they don't all 'die' on the
% same frame:
% DM.life = ceil(rand(1,DM.nDots)*DM.lifetime_fr);

for r = 1 : DM.n_reps
    switch Param.seqmode
        case 'random'
            DM.seqdirections(r,:) = randperm(DM.directions);
        case 'sequential'
            DM.seqdirections(r,:) = [1:DM.directions];  % to have sequential sequence of directions, not random.
    end
end

DM.seqangles = [DM.seqdirections * 360/DM.directions] + DM.offset_rot_deg ;  % random sequence of directions (e.g. random sequence of 30:30:360)
% list of angles given in cartesian coordinates (0deg=upward,
% increasing clockwise):
angles_cartesian = [1:DM.directions]*360/DM.directions;
DM.angles_cartesian = mod(angles_cartesian, 360);

Param.stimSeq(end+1,1) = { repmat(DM.stim_id, 1,DM.directions*DM.n_reps) };
DM.cnt=1;


if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
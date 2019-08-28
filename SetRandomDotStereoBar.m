function [ RDSB, Param ] = SetRandomDotStereoBar( RDSB, Param )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

%% settings for dots - taken from SetRandomDotStereogram

RDSB.apertureSize = []; % size of rectangular aperture [w,h] in degrees; if [], aperture=screen.
if isempty(RDSB.apertureSize)
    RDSB.apertureSize = fix([Param.screenRect(3),Param.screenRect(4)]/mean(Param.pixperdeg));
end
RDSB.center = [0,0]; % center of the field of dots (x,y)

RDSB.size_px = round( RDSB.size * Param.pixperdeg );
areaScreen = prod(Param.screenRes);
areaDot = (mean(RDSB.size_px)/2)^2*pi;
RDSB.nDots   = []; % number of dots; if empty nDots is calculated from the density.
if isempty(RDSB.nDots)
    RDSB.nDots = round( areaScreen*RDSB.density/areaDot );
end

% if     length(RDSB.color)==2 && ~isempty(RDSB.color{2})
%     %half nDots with color{1} and half with color{2}:
%     colorlist1 = repmat(RDSB.color{1}', 1,ceil(RDSB.nDots/2));
%     colorlist2 = repmat(RDSB.color{2}', 1, fix(RDSB.nDots/2));
%     RDSB.color_list{1} = uint8([colorlist1, colorlist2]);
%     if RDSB.Anticorrelated == 1
%     RDSB.color_list{2} = uint8([colorlist2, colorlist1]);
%     else
%     RDSB.color_list{2} = uint8([colorlist1, colorlist2]); 
%     end
% elseif length(RDSB.color)==1 || isempty(RDSB.color{2})
%     RDSB.color_list{1} = uint8(repmat(RDSB.color{1}', 1,RDSB.nDots));
%     %%% the nDots in the right screen is not exactly the same as in the
%     %%% left screen, so you need to get color_list{2} in PresentRDSB
%     %%% function
% %     RDSB.color_list{2} = RDSB.color_list{1};
% end

% We need to deal with DM moving beyond the edge of the aperture. This
% requrires a couple more lines of code.
% First we'll calculate the left,
% right top and bottom of the aperture (in degrees)
aperture.l = RDSB.center(1)-RDSB.apertureSize(1)/2;
aperture.r = RDSB.center(1)+RDSB.apertureSize(1)/2;
aperture.b = RDSB.center(2)-RDSB.apertureSize(2)/2;
aperture.t = RDSB.center(2)+RDSB.apertureSize(2)/2;
RDSB.aperture = aperture;

%%

if ~isempty(RDSB.angles) && numel(RDSB.angles)==RDSB.directions
    RDSB.seqangles = repmat(RDSB.angles-90, RDSB.n_disparities,1); % +90 to convert from cartesian to PTB angles
    angles_cartesian = RDSB.angles;
else
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    angles_cartesian = [1:RDSB.directions]*360/RDSB.directions;
    angles_cartesian = sort( mod(angles_cartesian, 360) );
    RDSB.seqangles = repmat(mod(angles_cartesian-90, 360), RDSB.n_disparities,1);
end
RDSB.seqangles(RDSB.seqangles == 360) = 0;
RDSB.angles_cartesian = mod(angles_cartesian, 360);
RDSB.angles_ptb = RDSB.angles_cartesian - 90;

RDSB.disparity_values_px = round( RDSB.disparity_values * Param.pixperdeg(1) );

RDSB.pattern_time_fr = round(RDSB.pattern_time/Param.ifi); % # frames per pattern.
RDSB.frames_poststim = round(RDSB.poststim_time/Param.ifi);
% to position the stimulation in the screen. Center of the screen:
RDSB.mouseCenter(1,:) = [ (Param.screenRect(1,3)-Param.screenRect(1,1))/2  (Param.screenRect(1,4)-Param.screenRect(1,2))/2];% in pixel coordinates
RDSB.mouse_dist_px(1) = fix((Param.screenRect(1,3) / Param.screenSize_cm_horz) * Param.mouse_dist_cm(1));

RDSB.bar_width_px = round( RDSB.bar_width  * Param.pixperdeg(1) );
RDSB.bar_lenth_px = round( RDSB.bar_length * Param.pixperdeg(1) );

RDSB.bar_speed_px = round( RDSB.bar_speed * Param.pixperdeg(1) ); %px/s
RDSB.bar_shift_perPattern = round( RDSB.bar_speed * RDSB.pattern_time );

% x = meshgrid(1:Param.screenRect(3), 1);
% RDSB.freq_cppx = RDSB.spacfreq/Param.pixperdeg(1);
% RDSB.phaseincrement = (RDSB.cyclespersecond * 360) * Param.ifi;

nTexBars =  round( RDSB.stimulus_time_perCycle/RDSB.pattern_time );
RDSB.nTexBars = nTexBars;

%convert to radians
bar_width_rad  = RDSB.bar_width / 180 * pi;
bar_length_rad = RDSB.bar_length / 180 * pi;
size_fov_rad   = RDSB.fov ./ 180 .* pi;
offset_fov_rad = (RDSB.view_offset - RDSB.fov/2) ./ 180 .* pi;
% we make only vertical bars (later we change the angle of the texture if we
% want orientation different than vertical):
n_stim = [nTexBars, 1];

%calculate the center x,y points of the patches
start_cent = offset_fov_rad + size_fov_rad./n_stim./2;  % center first patch
inc_cent = size_fov_rad./n_stim;                    % distance patch centers
end_cent = offset_fov_rad + size_fov_rad - size_fov_rad./n_stim./2;  % center last patch
cent_lo = start_cent(1):inc_cent(1):end_cent(1);  % centers along x axis
cent_la = start_cent(2):inc_cent(2):end_cent(2);  % centers along y axis

% vectors with all corner points:
long = zeros(prod(n_stim),4);
lat  = zeros(prod(n_stim),4);
i=1;
for x=1:n_stim(1) 
    for y=1:n_stim(2)
        long(i, [1,4]) = cent_lo(x) - bar_width_rad/2; 
        long(i, [2,3]) = cent_lo(x) + bar_width_rad/2; 
        lat(i, [1,2])  = cent_la(y) - bar_length_rad/2; 
        lat(i, [3,4])  = cent_la(y) + bar_length_rad/2; 
        i=i+1;
    end
end

% switch RDSB.gnomonic
%     case 1
%         [x,y]=  pr_gnomonic(reshape(long, [],1),reshape(lat, [],1));   %pr_gnomonic: coordinates with the Gnomonic non conformal projection; 
%         xx = reshape(x,[],4);
%         yy = reshape(y,[],4);
%     case 0
        xx(:,:) = long;
        yy(:,:) = lat;
% end

coord_x = round( xx.* RDSB.mouse_dist_px(1) + RDSB.mouseCenter(1,1) ); 
coord_y = round( yy.* RDSB.mouse_dist_px(1) + RDSB.mouseCenter(1,2) ); 


coord_x = coord_x(:,[1,2]);
coord_y = coord_y(:,[1,3]);
% % coord_x_left  = zeros(nTexBars, 2, RDSB.n_disparities, 'uint16');
% % coord_x_right = zeros(nTexBars, 2, RDSB.n_disparities, 'uint16');
% % for d = 1 : RDSB.n_disparities
% %     disp_shift_px = round( RDSB.disparity_values_px(d)/2 );
% %     coord_x_left(:,:,d)  = coord_x + disp_shift_px ;
% %     coord_x_right(:,:,d) = coord_x - disp_shift_px ;
% % end
% % coord_y_left  = repmat(coord_y, 1,1,RDSB.n_disparities);
% % coord_y_right = coord_y_left;



% extend screen by 2*marg to account for rotation:
if any(RDSB.angles_ptb ~= 0 | RDSB.angles_ptb ~= 180)
    marg = 1000;
else
    marg = 0;
end
screenExtended = gpuArray( false(Param.screenRect(1,4)+2*marg, Param.screenRect(1,3)+2*marg) );

coord_x_ext = coord_x + marg;
coord_y_ext = coord_y + marg;

mask_bar  = false(Param.screenRect(1,4), Param.screenRect(1,3), n_stim(1), 2, length(RDSB.angles_ptb), length(RDSB.disparity_values_px));

% figure;
for x = 1 : n_stim(1)
    screenExtended_tmp = (screenExtended);
    screenExtended_tmp([coord_y_ext(x,1):coord_y_ext(x,2)],[coord_x_ext(x,1):coord_x_ext(x,2)]) = true;
    for a = 1 : length(RDSB.angles_ptb)
        angle_ptb  = RDSB.angles_ptb(a);
        angle_cart = RDSB.angles_cartesian(a);
        mask_ext_rot = imrotate(screenExtended_tmp, -angle_ptb, 'bilinear','crop');
%         mask_rot     = mask_ext_rot(marg+1:end-marg,marg+1:end-marg);
        for d = 1 : length(RDSB.disparity_values_px)
%             if ismember(angle_ptb, [90, -90, 270]) % no gap with horizontal edge
%             else
                disparity_px = RDSB.disparity_values_px(d);
%                 if angle_cart > 180
%                     disphalf_px = -round( disparity_px/2 );
%                 else
                    disphalf_px = +round( disparity_px/2 );
%                 end
%             end
            
            mask_ext_left     = circshift( mask_ext_rot, [0 -disphalf_px]);
            mask_ext_right    = circshift( mask_ext_rot, [0 +disphalf_px]);
            mask_bar(:,:,x,1,a,d) = gather( mask_ext_left (marg+1:end-marg,marg+1:end-marg) );
            mask_bar(:,:,x,2,a,d) = gather( mask_ext_right(marg+1:end-marg,marg+1:end-marg) );
            
        end
    end
end

RDSB.mask_bar = mask_bar;

RDSB.coord_x = coord_x;
RDSB.coord_y = coord_y;
% RDSB.coord_x_left  = coord_x_left;
% RDSB.coord_x_right = coord_x_right;
% RDSB.coord_y_left  = coord_y_left;
% RDSB.coord_y_right = coord_y_right;


RDSB.seqdirections = repmat(1:RDSB.directions, RDSB.n_disparities,1);
RDSB.seqIntervals  = repmat([1:RDSB.n_disparities]', 1,RDSB.directions);
for r = 1 : RDSB.n_reps
    switch Param.seqmode
        case 'random'
            RDSB.stimseq(r,:) = randperm(RDSB.directions*RDSB.n_disparities);
        case 'sequential'
            RDSB.stimseq(r,:) = 1 : RDSB.directions*RDSB.n_disparities; 
    end
end

Param.stimSeq(end+1,1) = { repmat(RDSB.stim_id, 1,RDSB.directions*RDSB.n_disparities*RDSB.n_reps) };
RDSB.cnt=1;


if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end




function [ RDSB, Param ] = SetRandomDotStereoBar( RDSB, Param, win )


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

nTexBars =  round( RDSB.stimulus_time/RDSB.pattern_time );
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

switch RDSB.gnomonic
    case 1
        [x,y]=  pr_gnomonic(reshape(long, [],1),reshape(lat, [],1));   %pr_gnomonic: coordinates with the Gnomonic non conformal projection; 
        xx = reshape(x,[],4);
        yy = reshape(y,[],4);
    case 0
        xx(:,:) = long;
        yy(:,:) = lat;
end

coord_x = round( xx.* RDSB.mouse_dist_px(1) + RDSB.mouseCenter(1,1) ); 
coord_y = round( yy.* RDSB.mouse_dist_px(1) + RDSB.mouseCenter(1,2) ); 
coord_x = coord_x(:,[1,2]);
coord_y = coord_y(:,[1,3]);
coord_x_left  = zeros(nTexBars, 2, RDSB.n_disparities, 'uint16');
coord_x_right = zeros(nTexBars, 2, RDSB.n_disparities, 'uint16');
for d = 1 : RDSB.n_disparities
    disp_shift_px = round( RDSB.disparity_values_px(d)/2 );
    coord_x_left(:,:,d)  = coord_x + disp_shift_px ;
    coord_x_right(:,:,d) = coord_x - disp_shift_px ;
end
coord_y_left  = repmat(coord_y, 1,1,RDSB.n_disparities);
coord_y_right = coord_y_left;

RDSB.coord_x_left  = coord_x_left;
RDSB.coord_x_right = coord_x_right;
RDSB.coord_y_left  = coord_y_left;
RDSB.coord_y_right = coord_y_right;


% polymask = false( Param.screenRect(1,4), Param.screenRect(1,3), nTexBars);
% for t = 1 : nTexBars
% 
%     polymask(:,:,t) = poly2mask(xx(t,:).* RDSB.mouse_dist_px(1) + RDSB.mouseCenter(1,1), yy(t,:).* RDSB.mouse_dist_px(1) + RDSB.mouseCenter(1,2), Param.screenRect(1,4)-Param.screenRect(1,2), Param.screenRect(1,3)-Param.screenRect(1,1));
%     % polymask: 1 inside the bar, 0 outside.
%     % Multiplying by mouse_dist_px you take into account the distance from
%     % the screen and rescale all the measures in pixels to obtain the
%     % proper values in degrees.
% 
% end
% RDSB.polymask = polymask;

% tic
nrTotDots       = Param.screenRect(1,4) * Param.screenRect(1,3);
nrBrightDots = round( RDSB.pixel_density * nrTotDots );
% List_BrightPixels = zeros( nrBrightDots, 100, 'uint32');
% tot_patterns = 100;
% for i = 1 : tot_patterns
%     List_BrightPixels(:,i) =  randperm(nrTotDots, nrBrightDots)' ;
% end
% toc

% RDSB.List_BrightPixels = List_BrightPixels;

RDSB.Textures = zeros(nTexBars, RDSB.n_disparities, 2);

for d = 1 : RDSB.n_disparities
    
    for t = 1 : nTexBars

        
        dotScreen_left  = zeros( Param.screenRect(1,4), Param.screenRect(1,3), 'uint8' ) + RDSB.black;
        List_BrightPixels =  randperm(nrTotDots, nrBrightDots)';
        dotScreen_left(  List_BrightPixels ) = RDSB.white;
        dotScreen_left = reshape(dotScreen_left, Param.screenRect(1,4), Param.screenRect(1,3));
        dotScreen_right = zeros( Param.screenRect(1,4), Param.screenRect(1,3), 'uint8' ) + RDSB.black;
        List_BrightPixels =  randperm(nrTotDots, nrBrightDots)';
        dotScreen_right( List_BrightPixels  ) = RDSB.white;
        dotScreen_right = reshape(dotScreen_right, Param.screenRect(1,4), Param.screenRect(1,3));

        dotScreen_right([coord_y_right(t,1,d):coord_y_right(t,2,d)],[coord_x_right(t,1,d):coord_x_right(t,2,d)]) =...
            dotScreen_left([coord_y_left(t,1,d):coord_y_left(t,2,d)],[coord_x_left(t,1,d):coord_x_left(t,2,d)]);

        RDSB.Textures(t,d,1) = Screen('MakeTexture', win, dotScreen_left);
        RDSB.Textures(t,d,2) = Screen('MakeTexture', win, dotScreen_right);

    end
end


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

Param.stimSeq(end+1,1) = { repmat(RDSB.stim_id, 1,RDSB.directions*RDSB.n_disparities*RDSB.n_reps) };
RDSB.cnt=1;

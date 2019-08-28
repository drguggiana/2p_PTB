function [ RDSE, Param ] = SetRandomDotStereoEdge( RDSE, Param )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

%% settings for dots - taken from SetRandomDotStereogram

RDSE.apertureSize = []; % size of rectangular aperture [w,h] in degrees; if [], aperture=screen.
if isempty(RDSE.apertureSize)
    RDSE.apertureSize = fix([Param.screenRect(3),Param.screenRect(4)]/mean(Param.pixperdeg));
end
RDSE.center = [0,0]; % center of the field of dots (x,y)

RDSE.size_px = round( RDSE.size * Param.pixperdeg );
areaScreen = prod(Param.screenRes);
dotRadius = round(mean(RDSE.size_px)/2);
areaDot   = dotRadius^2*pi;
RDSE.nDots   = []; % number of dots; if empty nDots is calculated from the density.
if isempty(RDSE.nDots)
    RDSE.nDots = round( areaScreen*RDSE.density/areaDot );
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
aperture.l = RDSE.center(1)-RDSE.apertureSize(1)/2;
aperture.r = RDSE.center(1)+RDSE.apertureSize(1)/2;
aperture.b = RDSE.center(2)-RDSE.apertureSize(2)/2;
aperture.t = RDSE.center(2)+RDSE.apertureSize(2)/2;
RDSE.aperture = aperture;

%%

if ~isempty(RDSE.angles) && numel(RDSE.angles)==RDSE.directions
    RDSE.seqangles = repmat(RDSE.angles-90, RDSE.n_disparities,1); % +90 to convert from cartesian to PTB angles
    angles_cartesian = RDSE.angles;
else
    % list of angles given in cartesian coordinates (0deg=upward,
    % increasing clockwise):
    angles_cartesian = [1:RDSE.directions]*360/RDSE.directions;
    angles_cartesian = sort( mod(angles_cartesian, 360) );
    RDSE.seqangles = repmat(mod(angles_cartesian-90, 360), RDSE.n_disparities,1);
end
RDSE.seqangles(RDSE.seqangles == 360) = 0;
RDSE.angles_cartesian = mod(angles_cartesian, 360);
RDSE.angles_ptb = RDSE.angles_cartesian - 90;


RDSE.disparity_values_px = round( RDSE.disparity_values * Param.pixperdeg(1) );


RDSE.pattern_time_fr = round(RDSE.pattern_time/Param.ifi); % # frames per pattern.
RDSE.frames_poststim = round(RDSE.poststim_time/Param.ifi);
% to position the stimulation in the screen. Center of the screen:
RDSE.mouseCenter(1,:) = [ (Param.screenRect(1,3)-Param.screenRect(1,1))/2  (Param.screenRect(1,4)-Param.screenRect(1,2))/2];% in pixel coordinates
RDSE.mouse_dist_px(1) = fix((Param.screenRect(1,3) / Param.screenSize_cm_horz) * Param.mouse_dist_cm(1));

RDSE.bar_speed_px = round( RDSE.bar_speed * Param.pixperdeg(1) ); %px/s
RDSE.bar_shift_perPattern = round( RDSE.bar_speed * RDSE.pattern_time );

% x = meshgrid(1:Param.screenRect(3), 1);
% RDSB.freq_cppx = RDSB.spacfreq/Param.pixperdeg(1);
% RDSB.phaseincrement = (RDSB.cyclespersecond * 360) * Param.ifi;

nTexBars =  round( RDSE.stimulus_time_perCycle/RDSE.pattern_time );
RDSE.nTexBars = nTexBars;

%convert to radians
size_fov_rad   = RDSE.fov ./ 180 .* pi;
offset_fov_rad = (RDSE.view_offset - RDSE.fov/2) ./ 180 .* pi;
% we make only vertical bars (later we change the angle of the texture if we
% want orientation different than vertical):
n_stim = [nTexBars, 1];



%calculate the center x,y points of the patches
start_cent = offset_fov_rad + size_fov_rad./n_stim./2;  % center first patch
inc_cent   = size_fov_rad./n_stim;                    % distance patch centers
end_cent   = offset_fov_rad + size_fov_rad - size_fov_rad./n_stim./2;  % center last patch
cent_lo    = start_cent(1):inc_cent(1):end_cent(1);  % centers along x axis
cent_la    = start_cent(2):inc_cent(2):end_cent(2);  % centers along y axis

% vectors with all corner points:
long = zeros(prod(n_stim),4);
lat  = zeros(prod(n_stim),4);
i=1;
for x = 1 : n_stim(1) 
    for y=1:n_stim(2)
        long(i, 1)     = cent_lo(x); 
        lat (i, [1,2]) = [cent_la(y) - size_fov_rad(2)/2,...
                          cent_la(y) + size_fov_rad(2)/2]; 
        i=i+1;
    end
end

% switch RDSB.gnomonic
%     case 1
%         [x,y]=  pr_gnomonic(reshape(long, [],1),reshape(lat, [],1));   %pr_gnomonic: coordinates with the Gnomonic non conformal projection; 
%         xx = reshape(x,[],4);
%         yy = reshape(y,[],4);
%     case 0
        xx = long;
        yy = lat;
% end

coord_x = round( xx.* RDSE.mouse_dist_px(1) + RDSE.mouseCenter(1,1) ); 
coord_y = round( yy.* RDSE.mouse_dist_px(1) + RDSE.mouseCenter(1,2) ); 

% extend screen by 2*marg to account for rotation:
if any(RDSE.angles_ptb ~= 0 | RDSE.angles_ptb ~= 180)
    marg = 1000;
else
    marg = 0;
end
screenExtended = gpuArray( false(Param.screenRect(1,4)+2*marg, Param.screenRect(1,3)+2*marg) );

coord_x_ext = coord_x + marg;

mask_edge = false(Param.screenRect(1,4), Param.screenRect(1,3), n_stim(1), length(RDSE.angles_ptb));
mask_gap  = false(Param.screenRect(1,4), Param.screenRect(1,3), n_stim(1), length(RDSE.angles_ptb), length(RDSE.disparity_values_px));

% figure;
for x = 1 : n_stim(1)
    screenExtended_tmp = screenExtended;
    screenExtended_tmp(:,[1:coord_x_ext(x)]) = true;
%     screenExtended_le_tmp = screenExtended;
%     screenExtended_le_tmp(:,[1:coord_x_ext(x)-dotRadius]) = true;
%     screenExtended_re_tmp = screenExtended;
%     screenExtended_re_tmp(:,[1:coord_x_ext(x)+dotRadius]) = true;
%     screenExtended_re_tmp = ~screenExtended_re_tmp;
    for a = 1 : length(RDSE.angles_ptb)
        angle_ptb  = RDSE.angles_ptb(a);
        angle_cart = RDSE.angles_cartesian(a);
        mask_ext    = imrotate(screenExtended_tmp, -angle_ptb, 'bilinear','crop');
%         mask_ext_le = imrotate(screenExtended_le_tmp, -angle_ptb, 'bilinear','crop');
%         mask_ext_re = imrotate(screenExtended_re_tmp, -angle_ptb, 'bilinear','crop');
%         mask     = gather( mask_ext(marg+1:end-marg,marg+1:end-marg) );
%         mask_edge(:,:,x,1,a)     = gather( mask_ext_le(marg+1:end-marg,marg+1:end-marg) );
%         mask_edge(:,:,x,2,a)     = gather( mask_ext_re(marg+1:end-marg,marg+1:end-marg) );
        mask_edge(:,:,x,a)     = gather( mask_ext(marg+1:end-marg,marg+1:end-marg) );
%         imagesc(mask); drawnow;
%         [xy_le{x,a}(:,2), xy_le{x,a}(:,1)] = find(  mask );
%         [xy_re{x,a}(:,2), xy_re{x,a}(:,1)] = find( ~mask );
%         [ind_le{x,a}] = find(  mask );
%         [ind_re{x,a}] = find( ~mask );
%         im = false(Param.screenRes)'; iii = sub2ind(size(mask),y_le{x,a},x_le{x,a}); im(iii) = true; imagesc(im); drawnow;
        % find gap (it changes depending on the angle (e.g. horizontal edge
        % there is no gap)):
        for d = 1 : length(RDSE.disparity_values_px)
            if ismember(angle_ptb, [90, -90, 270]) % no gap with horizontal edge
%                 ind_gap{x,a,d} = [];
%                 mask_gap(:,:,x,a,d) = false(Param.screenRect(1,4), Param.screenRect(1,3));
            else
                disparity_px = RDSE.disparity_values_px(d);
%                 disp_px = round( disparity_px * cos(angle*pi/180) );
                if angle_cart > 180
                    disphalf_px = -round( disparity_px/2 );
                else
                    disphalf_px = +round( disparity_px/2 );
                end
                gap_ext_le = circshift( mask_ext, [0 -disphalf_px]);
                gap_ext_re = circshift(~mask_ext, [0 +disphalf_px]);
%                 gap_ext_le = circshift( mask_ext_le, [0 -disphalf_px]);
%                 gap_ext_re = circshift( mask_ext_re, [0 +disphalf_px]);

                gap_le = gap_ext_le(marg+1:end-marg,marg+1:end-marg);
                gap_re = gap_ext_re(marg+1:end-marg,marg+1:end-marg);
                gap_sumlr = gather( gap_le + gap_re );
                gap_sumlr = logical(gap_sumlr);
                mask_gap(:,:,x,a,d) = ~gap_sumlr;
%                 imagesc(gap_sumlr); drawnow;
                % gap consists of false values in gap_sumlr:
%                 [xy_gap{x,a,d}(:,2), xy_gap{x,a,d}(:,1)] = find( ~gap_sumlr );
%                 [ind_gap{x,a,d}] = find( ~gap_sumlr );
%               im = false(Param.screenRes)'; iii = sub2ind(size(maskg), y_gap{x,a,d}, x_gap{x,a,d}); im(iii) = true; imagesc(im); drawnow;
            end
        end
    end
end

RDSE.mask_edge = mask_edge;
RDSE.mask_gap  = mask_gap;



RDSE.seqdirections = repmat(1:RDSE.directions, RDSE.n_disparities,1);
RDSE.seqIntervals  = repmat([1:RDSE.n_disparities]', 1,RDSE.directions);
for r = 1 : RDSE.n_reps
    switch Param.seqmode
        case 'random'
            RDSE.stimseq(r,:) = randperm(RDSE.directions*RDSE.n_disparities);
        case 'sequential'
            RDSE.stimseq(r,:) = 1 : RDSE.directions*RDSE.n_disparities; 
    end
end


Param.stimSeq(end+1,1) = { repmat(RDSE.stim_id, 1,RDSE.directions*RDSE.n_disparities*RDSE.n_reps) };
RDSE.cnt=1;


if exist('scriptName','var')
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end




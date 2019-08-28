[scriptPathPhi, scriptNamePhi] = fileparts(mfilename('fullpath'));

% Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% show dark screen while all the bars are generated...
Screen('FillRect', screen(1), 70);
Screen('Flip', screen(1));
if Param.Dichoptic
    Screen('FillRect', screen(2), 70);
    Screen('Flip', screen(2));
end

PM.positions = [1, 1];  % in how many x and y positions divide the FOV.
PM.rotation_rad = PM.rotation_deg * pi/180 ; % radians
PM.offset_rot_rad = PM.offset_rot_deg * pi/180 ; % radians
% PM.locations = length(PM.bar_centers) ;
PM.configurations = PM.locations * numel(PM.rotation_deg) ; % !! don't forget to consider black and white (not considered here)!
% PM.nbr_stimuli = (PM.Black+PM.White)*(PM.locations^2)*numel(PM.rotation_deg) ;
% if Param.Dichoptic
%     PM.nbr_stimuli = 2*PM.nbr_stimuli ; % there are additionalStaticBars
% end
PM.mouseCenter(1,:) = [ (Param.screenRect(1,3)-Param.screenRect(1,1))/2  (Param.screenRect(1,4)-Param.screenRect(1,2))/2];%19/25 ] ; % in pixel coordinates (position the mouse pointer on the screen an use GetMouse in MatLab)
PM.mouse_dist_px(1) = fix((Param.screenRect(1,3) / screenSize_cm_horz) * mouse_dist_cm(1));% * 0.8);  % in pixel
    % I added the factor 0.8 because for some reason it does not scale properly the bars (e.g. 12 deg is actually presented as 14.8) 
if Param.Dichoptic
    PM.mouseCenter(2,:) = [ (Param.screenRect(1,3)-Param.screenRect(1,1))/2  (Param.screenRect(1,4)-Param.screenRect(1,2))/2];
    PM.mouse_dist_px(2) = fix((Param.screenRect(1,3) / screenSize_cm_horz) * mouse_dist_cm(2));% * 0.8);  % in pixel
end
PM.frame_stim = round(PM.patch_time/Param.ifi);          %Ale: # frames per patch.
% PM.frame_poststim = round(PM.poststim_time/Param.ifi);  %Ale: # frames per poststim blank screen.
PM.frame_interpatch = round(PM.interpatch_time/Param.ifi);

% to smooth the edge
if PM.edgeSmoothing
    widthBar_px = PM.widthBar_deg * Param.pixperdeg(1);
    fraction = 4;
    hsize = [round( widthBar_px/fraction), round(widthBar_px/fraction)]
    % sigma = 6;
    % h = fspecial('gaussian', hsize, sigma);
    h = fspecial('average', hsize);
    % the smoothing enlarges the bar, so bar width needs to be reduced
    % accordingly:
    widthBar_deg = PM.widthBar_deg /(1+1/fraction);
end

n = 1;
for xx = 1 : PM.positions(1)
position_L_x_cent(xx) = n * PM.fov(1)/PM.positions(1)/2 + PM.offset_stim(1,1);
if Param.Dichoptic
    position_R_x_cent(xx) = n * PM.fov(1)/PM.positions(1)/2 + PM.offset_stim(2,1);
else
    position_R_x_cent(xx) = position_L_x_cent(xx) ;
end
n = n+2;
end
n = 1;
for yy = 1 : PM.positions(2)
position_L_y_cent(yy) = n * PM.fov(2)/PM.positions(2)/2 + PM.offset_stim(1,2);
if Param.Dichoptic
    position_R_y_cent(yy) = n * PM.fov(2)/PM.positions(2)/2 + PM.offset_stim(2,2);
else
    position_R_y_cent(yy) = position_L_y_cent(yy) ;
end
n = n+2;
end

%%% pre-allocation of variables %%%
PM.masktex   = zeros(  PM.configurations  ,  2 );  % 2=2 bars, 12=configurations: 2x3x2 (2locations x 3distances x 2directions).
PM.polymask  = cell(   PM.configurations ,  2  , PM.positions(1) );
if PM.save_screen==1
    save_dir_maskscreen = [save_dir filesep 'screen_masks_' datetime_suffix filesep];
    if ~isdir(save_dir_maskscreen)
        mkdir(save_dir_maskscreen);
    end
end
if exist('RetinotopicPosition','var') && RetinotopicPosition
    patches=RP.tot_positions;
else
    patches=[10,8];
end
for i=1:length(screen)
    ScreenImg_merge(:,:,:,i) = DrawScreen(Param.screenRect,patches);
end
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

for p = 1:PM.positions(1)  %from position 1 to 4 calculate the coordinates of all bars.

%%% pre-allocation of variables for speed %%%
x_L_cent    = zeros( PM.locations , 1) ; % 4=number of distances+zero_distance, 1=only one coordinate (x of the centre of the first bar).
y_L_cent    = x_L_cent ; %
x_R_cent    = x_L_cent ; %
y_R_cent    = x_L_cent ; %
x_L_corners = zeros( PM.locations , 4 ) ;
y_L_corners = x_L_corners ; %
x_R_corners = x_L_corners ; %
y_R_corners = x_L_corners ; %
xx          = zeros( PM.configurations , 4 ) ;
xx(:,:,2)   = xx ;
yy          = xx ;
% toclist     = zeros( 2*PM.configurations^2, 4, PM.n_reps ) ; % to monitor the timings. 4= 1 is prestim, 2and3 are two flashes, 4 is poststim.
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

 position = [ p , 1 ] ;

% xy coordinates of the center of the first bar
for ii = 1 : PM.locations
    x_L_cent(ii,1) = position_L_x_cent(position(1)) + PM.bar_centers(ii) ;  % centers at -15 -10 -5 0 +5 +10 +15
    x_R_cent(ii,1) = position_R_x_cent(position(1)) + PM.bar_centers(ii) ;
end
    y_L_cent(1:PM.locations,1) = position_L_y_cent(position(2)) ;  % deg, position in the field of view.
    y_R_cent(1:PM.locations,1) = position_R_y_cent(position(2)) ;

% produces vectors with corner points for the first bar
x_L_corners(:,:)  = [ x_L_cent-(widthBar_deg/2)  , x_L_cent+(widthBar_deg/2) , x_L_cent+(widthBar_deg/2)  , x_L_cent-(widthBar_deg/2)  ];
y_L_corners(:,:)  = [ y_L_cent+(PM.lengthBar_deg/2) , y_L_cent+(PM.lengthBar_deg/2), y_L_cent-(PM.lengthBar_deg/2) , y_L_cent-(PM.lengthBar_deg/2) ];
x_R_corners(:,:)  = [ x_R_cent-(widthBar_deg/2)  , x_R_cent+(widthBar_deg/2) , x_R_cent+(widthBar_deg/2)  , x_R_cent-(widthBar_deg/2)  ];
y_R_corners(:,:)  = [ y_R_cent+(PM.lengthBar_deg/2) , y_R_cent+(PM.lengthBar_deg/2), y_R_cent-(PM.lengthBar_deg/2) , y_R_cent-(PM.lengthBar_deg/2) ];


%%% rotation of the bars %%%

% translate x and y to set x_cent and y_cent as origin and then apply
% rotation:
x_L_trans  = x_L_corners - position_L_x_cent(position(1)) ; % x  = [ x(1)-x_cent , x(2)-x_cent , x(3)-x_cent , x(4)-x_cent ];
y_L_trans  = y_L_corners - position_L_y_cent(position(2)) ; % y  = [ y(1)-y_cent , y(2)-y_cent , y(3)-y_cent , y(4)-y_cent ];
x_R_trans  = x_R_corners - position_R_x_cent(position(1)) ; % x  = [ x(1)-x_cent , x(2)-x_cent , x(3)-x_cent , x(4)-x_cent ];
y_R_trans  = y_R_corners - position_R_y_cent(position(2)) ; % y  = [ y(1)-y_cent , y(2)-y_cent , y(3)-y_cent , y(4)-y_cent ];

% apply rotation
for r = 1:numel(PM.rotation_deg)  % from 1 to 2 (number of directions).
% for d = 1:(distances-1) % for the 3 distances between the two flashed bars.
x_L_rot(:,:,r) = (x_L_trans(:,:) * cos(PM.rotation_rad(r)+PM.offset_rot_rad(1))) - (y_L_trans(:,:) * sin(PM.rotation_rad(r)+PM.offset_rot_rad(1))) ;
y_L_rot(:,:,r) = (x_L_trans(:,:) * sin(PM.rotation_rad(r)+PM.offset_rot_rad(1))) + (y_L_trans(:,:) * cos(PM.rotation_rad(r)+PM.offset_rot_rad(1))) ;
x_R_rot(:,:,r) = (x_R_trans(:,:) * cos(PM.rotation_rad(r)+PM.offset_rot_rad(2))) - (y_R_trans(:,:) * sin(PM.rotation_rad(r)+PM.offset_rot_rad(2))) ;
y_R_rot(:,:,r) = (x_R_trans(:,:) * sin(PM.rotation_rad(r)+PM.offset_rot_rad(2))) + (y_R_trans(:,:) * cos(PM.rotation_rad(r)+PM.offset_rot_rad(2))) ;
end, % end
x_L_rot = permute(x_L_rot,[1 3 2]) ;
y_L_rot = permute(y_L_rot,[1 3 2]) ;
x_R_rot = permute(x_R_rot,[1 3 2]) ;
y_R_rot = permute(y_R_rot,[1 3 2]) ;
x_L_rot = reshape( x_L_rot , PM.configurations , 4 ) ;
y_L_rot = reshape( y_L_rot , PM.configurations , 4 ) ;
x_R_rot = reshape( x_R_rot , PM.configurations , 4 ) ;
y_R_rot = reshape( y_R_rot , PM.configurations , 4 ) ;

% reposition the origin of reference frame:
x_L_rot = x_L_rot + position_L_x_cent(position(1)) ;
y_L_rot = y_L_rot + position_L_y_cent(position(2)) ;
x_R_rot = x_R_rot + position_R_x_cent(position(1)) ;
y_R_rot = y_R_rot + position_R_y_cent(position(2)) ;

%%% rotation completed %%%

% to convert the coordinates you see in the RF outlines into correspondent
% coordinates of Screen('MakeTexture'):
lo_L =  (x_L_rot - PM.fov(1)/2) * pi/180 ;
la_L =  (y_L_rot - PM.fov(2)/2) * pi/180 ;
lo_R =  (x_R_rot - PM.fov(1)/2) * pi/180 ;
la_R =  (y_R_rot - PM.fov(2)/2) * pi/180 ;
% pr_gnomonic: coordinates with the Gnomonic non conformal projection:
switch PM.gnomonic
    case 1
        [x_L,y_L] = pr_gnomonic(reshape(lo_L,[],1), reshape(la_L,[],1)); % reshape lo and la in one single column
        [x_R,y_R] = pr_gnomonic(reshape(lo_R,[],1), reshape(la_R,[],1));      
        xx(:,:,1) = reshape(x_L, PM.configurations ,4 );   % restore the matrix with 4 columns
        yy(:,:,1) = reshape(y_L, PM.configurations ,4 );
        xx(:,:,2) = reshape(x_R, PM.configurations ,4 );
        yy(:,:,2) = reshape(y_R, PM.configurations ,4 );
    case 0
        xx(:,:,1) = lo_L;
        yy(:,:,1) = la_L;
        xx(:,:,2) = lo_R;
        yy(:,:,2) = la_R;
end


for f = 1:2  % the two flashed bars for each configuration.
for c = 1:PM.configurations % 
    % polymask: 1 inside the bar, 0 outside.
    % Multiplying by mouse_dist_px you take into account the distance from
    % the screen and rescale all the measures in pixels to obtain the
    % proper values in degrees.
    if Param.Dichoptic
        PM.polymask{c,f,p} = poly2mask(reshape(xx(c,:,f),1,[]) .* PM.mouse_dist_px(f) + PM.mouseCenter(f,1) , reshape(yy(c,:,f),1,[]) .* PM.mouse_dist_px(f) + PM.mouseCenter(f,2),Param.screenRect(1,4),Param.screenRect(1,3)-Param.screenRect(1,1)) ;
    else
        PM.polymask{c,f,p} = poly2mask(reshape(xx(c,:,f),1,[]) .* PM.mouse_dist_px(1) + PM.mouseCenter(1,1) , reshape(yy(c,:,f),1,[]) .* PM.mouse_dist_px(1) + PM.mouseCenter(1,2),Param.screenRect(1,4),Param.screenRect(1,3)-Param.screenRect(1,1)) ;
    end
    % BW are LA planes where A is alpha, the transparency of a pixel.
    % Alpha values range between zero(=fully transparent) and 255(=fully
    % opaque).
    BW       = zeros( Param.screenRect(1,4), Param.screenRect(1,3)-Param.screenRect(1,1) , 2 ) ;
    BW(:,:,1)= PM.BackgroundLuminance ; 
    % producing BW for bars in either position.  0 inside the patch, 255 outside.
    BW(:,:,2) = 255-255*PM.polymask{c,f,p} ;
    if PM.edgeSmoothing
        BWtemp = padarray(BW(:,:,2),[100 100], 255, 'both');
        BWtemp = imfilter(BWtemp, h);
        BW(:,:,2) = BWtemp(101:end-100, 101:end-100);
    end
    if Param.Dichoptic
        PM.masktex(c,f) = Screen('MakeTexture', screen(f), BW ) ;
    else
        PM.masktex(c,f) = Screen('MakeTexture', screen(1), BW ) ;
    end
%     if PM.save_screen==1
%         ScreenImg = DrawScreen(screenRect,[8,6]); %screenRect=[0,0,1600,900]
%         ScreenImg(:,:,1) = ScreenImg(:,:,1) - PM.polymask{c,f,p} ;
%         ScreenImg(:,:,3) = ScreenImg(:,:,1);
%         imwrite(ScreenImg,[save_dir_maskscreen 'screenmask_config' num2str(c) '.png'],'png');
%     end
    if Param.Dichoptic
            ScreenImg_merge(:,:,2,f) = ScreenImg_merge(:,:,2,f) - PM.polymask{c,f,p} ;
            ScreenImg_merge(:,:,3,f) = ScreenImg_merge(:,:,2,f);
    else
        ScreenImg_merge(:,:,2,1) = ScreenImg_merge(:,:,2,1) - PM.polymask{c,f,p} ;
        ScreenImg_merge(:,:,3,1) = ScreenImg_merge(:,:,2,1);
    end
%     end

end
end

end

if PM.save_screen==1
    for i=1:length(screen)
        imwrite(ScreenImg_merge(:,:,:,i),[save_dir_maskscreen 'screenmask' num2str(i) '_merge.png'],'png');
    end
end
% if show_screen_only==1
%     imshow(ScreenImg_merge);
%     return
% end

PM.seqbars=[];
[A,B]=meshgrid( 1:PM.locations , 1:PM.locations );
seqbarstmp = reshape( cat(2,A,B) ,[],2 );
if Param.Dichoptic
    additionalStaticBars = repmat([1,1; 2,2; 3,3] ,3,1);
    seqbarstmp=cat(1,seqbarstmp,additionalStaticBars);
end
for ii = 1:numel(PM.rotation_deg)
    PM.seqbars = [PM.seqbars ; seqbarstmp+(ii-1)*PM.locations ] ;
end
if PM.Black+PM.White == 2
    PM.seqbars = [PM.seqbars ; PM.seqbars ] ; % duplicate for black and white.
end
if     PM.Black == 1 && PM.White == 1
    PM.seqcolor( 1:PM.nbr_stimuli/2 ,1) = PM.BlackLuminance ;       % first half black
    PM.seqcolor( PM.nbr_stimuli/2 +1 : PM.nbr_stimuli ,1) = PM.WhiteLuminance ;   % second half white
elseif PM.Black == 1 && PM.White == 0
    PM.seqcolor( 1:PM.nbr_stimuli ,1) = PM.BlackLuminance ;
elseif PM.Black == 0 && PM.White == 1
    PM.seqcolor( 1:PM.nbr_stimuli ,1) = PM.WhiteLuminance ;
end
if Param.Dichoptic
%     PM.screensequence = [1 2; 1 2; 1 2; 2 1; 1 2; 1 2; 2 1; 2 1; 1 2;...
%                             1 1; 1 1; 1 1; 2 2; 2 2; 2 2; 2 1; 2 1; 2 1];
    PM.screensequence = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2;...
                            1 1; 1 1; 1 1; 2 2; 2 2; 2 2; 1 2; 1 2; 1 2];
    PM.screensequence = repmat(PM.screensequence, (PM.Black+PM.White)*numel(PM.rotation_deg),1);
else
    PM.screensequence = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
    PM.screensequence = repmat(PM.screensequence, (PM.Black+PM.White)*numel(PM.rotation_deg),1);
end
for r = 1:PM.n_reps
    % seqpositionstmp is [1 1; 1 2; 1 3; 1 4; 2 1; ...; 6 4]:
    switch Param.seqmode
        case 'sequential'
            PM.stimseq(r,:) = 1 : PM.nbr_stimuli ;
        case 'random'
            PM.stimseq(r,:) = randperm( PM.nbr_stimuli ) ;
    end     
end

Param.stimSeq(end+1,1) = { repmat(PM.stim_id, 1,PM.nbr_stimuli*PM.n_reps) };
% stimSeqtmp4 = [ repmat(PM.stim_id, 1,PM.nbr_stimuli*PM.n_reps) ];
PM.cnt=1;  

if PM.save_screen==1
    for s=1:PM.nbr_stimuli/(PM.Black+PM.White)
        ScreenImg = DrawScreen(Param.screenRect,patches); %screenRect=[0,0,1600,900]
        for f = [1,2]%PM.screensequence % [1,2]
            if     f==1%PM.screensequence(1) % first bar in light green
            ScreenImg(:,:,1) = ScreenImg(:,:,1) - PM.polymask{PM.seqbars(s,f),f,1} ;
            elseif f==2%PM.screensequence(2) % second bar in dark green
            ScreenImg(:,:,1) = ScreenImg(:,:,1) - 0.7*PM.polymask{PM.seqbars(s,f),f,1} ;
            ScreenImg(:,:,2) = ScreenImg(:,:,2) - 0.4*PM.polymask{PM.seqbars(s,f),f,1} ;
            end
            ScreenImg(:,:,3) = ScreenImg(:,:,1);
            imwrite(ScreenImg,[save_dir_maskscreen 'screenmask_config' num2str(s) '.png'],'png');
        end
    end
end

if exist('scriptNamePhi','var')
    BackupStimulusScript( scriptPathPhi, scriptNamePhi, save_dir )
end
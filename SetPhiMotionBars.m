[scriptPathPhi, scriptNamePhi] = fileparts(mfilename('fullpath'));

% Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% show dark screen while all the bars are generated...
% Screen('FillRect', screen(1), 70);
% Screen('Flip', screen(1));
% if Param.Dichoptic
%     Screen('FillRect', screen(2), 70);
%     Screen('Flip', screen(2));
% end
Screen('FillRect', screen, 70);

PMB.positions = [1, 1];  % in how many x and y positions divide the FOV.
PMB.rotation_rad = PMB.rotation_deg * pi/180 ; % radians
PMB.offset_rot_rad = PMB.offset_rot_deg * pi/180 ; % radians
% PMB.locations = length(PMB.bar_centers) ;
PMB.configurations = PMB.locations * numel(PMB.rotation_deg) ; % !! don't forget to consider black and white (not considered here)!
% PMB.nbr_stimuli = (PMB.Black+PMB.White)*(PMB.locations^2)*numel(PMB.rotation_deg) ;
% if Param.Dichoptic
%     PMB.nbr_stimuli = 2*PMB.nbr_stimuli ; % there are additionalStaticBars
% end
PMB.mouseCenter(1,:) = [ (Param.screenRect(1,3)-Param.screenRect(1,1))/2  (Param.screenRect(1,4)-Param.screenRect(1,2))/2];%19/25 ] ; % in pixel coordinates (position the mouse pointer on the screen an use GetMouse in MatLab)
PMB.mouse_dist_px(1) = fix((Param.screenRect(1,3) / screenSize_cm_horz) * mouse_dist_cm(1));% * 0.8);  % in pixel
    % I added the factor 0.8 because for some reason it does not scale properly the bars (e.g. 12 deg is actually presented as 14.8) 
if Param.Dichoptic
    PMB.mouseCenter(2,:) = [ (Param.screenRect(1,3)-Param.screenRect(1,1))/2  (Param.screenRect(1,4)-Param.screenRect(1,2))/2];
    PMB.mouse_dist_px(2) = fix((Param.screenRect(1,3) / screenSize_cm_horz) * mouse_dist_cm(2));% * 0.8);  % in pixel
end
PMB.frame_stim = round(PMB.patch_time/Param.ifi);          %Ale: # frames per patch.
% PMB.frame_poststim = round(PMB.poststim_time/Param.ifi);  %Ale: # frames per poststim blank screen.
PMB.frame_interpatch = round(PMB.interpatch_time/Param.ifi);

% to smooth the edge
if PMB.edgeSmoothing
    widthBar_px = PMB.widthBar_deg * Param.pixperdeg(1);
    fraction = 4;
    hsize = [round( widthBar_px/fraction), round(widthBar_px/fraction)]
    % sigma = 6;
    % h = fspecial('gaussian', hsize, sigma);
    h = fspecial('average', hsize);
    % the smoothing enlarges the bar, so bar width needs to be reduced
    % accordingly:
    widthBar_deg = PMB.widthBar_deg /(1+1/fraction);
else
    widthBar_deg = PMB.widthBar_deg;
end

n = 1;
for xx = 1 : PMB.positions(1)
position_L_x_cent(xx) = n * PMB.fov(1)/PMB.positions(1)/2 + PMB.offset_stim(1,1);
if Param.Dichoptic
    position_R_x_cent(xx) = n * PMB.fov(1)/PMB.positions(1)/2 + PMB.offset_stim(2,1);
else
    position_R_x_cent(xx) = position_L_x_cent(xx) ;
end
n = n+2;
end
n = 1;
for yy = 1 : PMB.positions(2)
position_L_y_cent(yy) = n * PMB.fov(2)/PMB.positions(2)/2 + PMB.offset_stim(1,2);
if Param.Dichoptic
    position_R_y_cent(yy) = n * PMB.fov(2)/PMB.positions(2)/2 + PMB.offset_stim(2,2);
else
    position_R_y_cent(yy) = position_L_y_cent(yy) ;
end
n = n+2;
end

%%% pre-allocation of variables %%%
PMB.masktex   = zeros(  PMB.configurations  ,  2 );  % 2=2 bars, 12=configurations: 2x3x2 (2locations x 3distances x 2directions).
PMB.polymask  = cell(   PMB.configurations ,  2  , PMB.positions(1) );
if PMB.save_screen==1
    save_dir_maskscreen = [save_dir filesep 'screen_masks_' datetime_suffix filesep];
    if ~isdir(save_dir_maskscreen)
        mkdir(save_dir_maskscreen);
    end
end
if exist('RetinotopicPosition','var') && RetinotopicPosition
    patches=RP.tot_positions;
else
    patches=[18,6];
end

ScreenImg_merge(:,:,:,1) = DrawScreen(Param.screenRect,patches);
ScreenImg_merge(:,:,:,2) = ScreenImg_merge(:,:,:,1);

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

for p = 1:PMB.positions(1)  %from position 1 to 4 calculate the coordinates of all bars.

%%% pre-allocation of variables for speed %%%
x_L_cent    = zeros( PMB.locations , 1) ; % 4=number of distances+zero_distance, 1=only one coordinate (x of the centre of the first bar).
y_L_cent    = x_L_cent ; %
x_R_cent    = x_L_cent ; %
y_R_cent    = x_L_cent ; %
x_L_corners = zeros( PMB.locations , 4 ) ;
y_L_corners = x_L_corners ; %
x_R_corners = x_L_corners ; %
y_R_corners = x_L_corners ; %
xx          = zeros( PMB.configurations , 4 ) ;
xx(:,:,2)   = xx ;
yy          = xx ;
% toclist     = zeros( 2*PMB.configurations^2, 4, PMB.n_reps ) ; % to monitor the timings. 4= 1 is prestim, 2and3 are two flashes, 4 is poststim.
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

 position = [ p , 1 ] ;

% xy coordinates of the center of the first bar
for ii = 1 : PMB.locations
    x_L_cent(ii,1) = position_L_x_cent(position(1)) + PMB.bar_centers(ii) ;  % centers at -15 -10 -5 0 +5 +10 +15
    x_R_cent(ii,1) = position_R_x_cent(position(1)) + PMB.bar_centers(ii) ;
end
    y_L_cent(1:PMB.locations,1) = position_L_y_cent(position(2)) ;  % deg, position in the field of view.
    y_R_cent(1:PMB.locations,1) = position_R_y_cent(position(2)) ;

% produces vectors with corner points for the first bar
x_L_corners(:,:)  = [ x_L_cent-(widthBar_deg/2)  , x_L_cent+(widthBar_deg/2) , x_L_cent+(widthBar_deg/2)  , x_L_cent-(widthBar_deg/2)  ];
y_L_corners(:,:)  = [ y_L_cent+(PMB.lengthBar_deg/2) , y_L_cent+(PMB.lengthBar_deg/2), y_L_cent-(PMB.lengthBar_deg/2) , y_L_cent-(PMB.lengthBar_deg/2) ];
x_R_corners(:,:)  = [ x_R_cent-(widthBar_deg/2)  , x_R_cent+(widthBar_deg/2) , x_R_cent+(widthBar_deg/2)  , x_R_cent-(widthBar_deg/2)  ];
y_R_corners(:,:)  = [ y_R_cent+(PMB.lengthBar_deg/2) , y_R_cent+(PMB.lengthBar_deg/2), y_R_cent-(PMB.lengthBar_deg/2) , y_R_cent-(PMB.lengthBar_deg/2) ];


%%% rotation of the bars %%%

% translate x and y to set x_cent and y_cent as origin and then apply
% rotation:
x_L_trans  = x_L_corners - position_L_x_cent(position(1)) ; % x  = [ x(1)-x_cent , x(2)-x_cent , x(3)-x_cent , x(4)-x_cent ];
y_L_trans  = y_L_corners - position_L_y_cent(position(2)) ; % y  = [ y(1)-y_cent , y(2)-y_cent , y(3)-y_cent , y(4)-y_cent ];
x_R_trans  = x_R_corners - position_R_x_cent(position(1)) ; % x  = [ x(1)-x_cent , x(2)-x_cent , x(3)-x_cent , x(4)-x_cent ];
y_R_trans  = y_R_corners - position_R_y_cent(position(2)) ; % y  = [ y(1)-y_cent , y(2)-y_cent , y(3)-y_cent , y(4)-y_cent ];

% apply rotation
for r = 1:numel(PMB.rotation_deg)  % from 1 to 2 (number of directions).
% for d = 1:(distances-1) % for the 3 distances between the two flashed bars.
x_L_rot(:,:,r) = (x_L_trans(:,:) * cos(PMB.rotation_rad(r)+PMB.offset_rot_rad(1))) - (y_L_trans(:,:) * sin(PMB.rotation_rad(r)+PMB.offset_rot_rad(1))) ;
y_L_rot(:,:,r) = (x_L_trans(:,:) * sin(PMB.rotation_rad(r)+PMB.offset_rot_rad(1))) + (y_L_trans(:,:) * cos(PMB.rotation_rad(r)+PMB.offset_rot_rad(1))) ;
x_R_rot(:,:,r) = (x_R_trans(:,:) * cos(PMB.rotation_rad(r)+PMB.offset_rot_rad(2))) - (y_R_trans(:,:) * sin(PMB.rotation_rad(r)+PMB.offset_rot_rad(2))) ;
y_R_rot(:,:,r) = (x_R_trans(:,:) * sin(PMB.rotation_rad(r)+PMB.offset_rot_rad(2))) + (y_R_trans(:,:) * cos(PMB.rotation_rad(r)+PMB.offset_rot_rad(2))) ;
end, % end
x_L_rot = permute(x_L_rot,[1 3 2]) ;
y_L_rot = permute(y_L_rot,[1 3 2]) ;
x_R_rot = permute(x_R_rot,[1 3 2]) ;
y_R_rot = permute(y_R_rot,[1 3 2]) ;
x_L_rot = reshape( x_L_rot , PMB.configurations , 4 ) ;
y_L_rot = reshape( y_L_rot , PMB.configurations , 4 ) ;
x_R_rot = reshape( x_R_rot , PMB.configurations , 4 ) ;
y_R_rot = reshape( y_R_rot , PMB.configurations , 4 ) ;

% reposition the origin of reference frame:
x_L_rot = x_L_rot + position_L_x_cent(position(1)) ;
y_L_rot = y_L_rot + position_L_y_cent(position(2)) ;
x_R_rot = x_R_rot + position_R_x_cent(position(1)) ;
y_R_rot = y_R_rot + position_R_y_cent(position(2)) ;

%%% rotation completed %%%

% to convert the coordinates you see in the RF outlines into correspondent
% coordinates of Screen('MakeTexture'):
lo_L =  (x_L_rot - PMB.fov(1)/2) * pi/180 ;
la_L =  (y_L_rot - PMB.fov(2)/2) * pi/180 ;
lo_R =  (x_R_rot - PMB.fov(1)/2) * pi/180 ;
la_R =  (y_R_rot - PMB.fov(2)/2) * pi/180 ;
% pr_gnomonic: coordinates with the Gnomonic non conformal projection:
switch PMB.gnomonic
    case 1
        [x_L,y_L] = pr_gnomonic(reshape(lo_L,[],1), reshape(la_L,[],1)); % reshape lo and la in one single column
        [x_R,y_R] = pr_gnomonic(reshape(lo_R,[],1), reshape(la_R,[],1));      
        xx(:,:,1) = reshape(x_L, PMB.configurations ,4 );   % restore the matrix with 4 columns
        yy(:,:,1) = reshape(y_L, PMB.configurations ,4 );
        xx(:,:,2) = reshape(x_R, PMB.configurations ,4 );
        yy(:,:,2) = reshape(y_R, PMB.configurations ,4 );
    case 0
        xx(:,:,1) = lo_L;
        yy(:,:,1) = la_L;
        xx(:,:,2) = lo_R;
        yy(:,:,2) = la_R;
end


for f = 1:2  % the two flashed bars for each configuration.
for c = 1:PMB.configurations % 
    % polymask: 1 inside the bar, 0 outside.
    % Multiplying by mouse_dist_px you take into account the distance from
    % the screen and rescale all the measures in pixels to obtain the
    % proper values in degrees.
    if Param.Dichoptic
        PMB.polymask{c,f,p} = poly2mask(reshape(xx(c,:,f),1,[]) .* PMB.mouse_dist_px(f) + PMB.mouseCenter(f,1) , reshape(yy(c,:,f),1,[]) .* PMB.mouse_dist_px(f) + PMB.mouseCenter(f,2),Param.screenRect(1,4),Param.screenRect(1,3)-Param.screenRect(1,1)) ;
    else
        PMB.polymask{c,f,p} = poly2mask(reshape(xx(c,:,f),1,[]) .* PMB.mouse_dist_px(1) + PMB.mouseCenter(1,1) , reshape(yy(c,:,f),1,[]) .* PMB.mouse_dist_px(1) + PMB.mouseCenter(1,2),Param.screenRect(1,4),Param.screenRect(1,3)-Param.screenRect(1,1)) ;
    end
    % BW are LA planes where A is alpha, the transparency of a pixel.
    % Alpha values range between zero(=fully transparent) and 255(=fully
    % opaque).
    BW       = zeros( Param.screenRect(1,4), Param.screenRect(1,3)-Param.screenRect(1,1) , 2 ) ;
    BW(:,:,1)= PMB.BackgroundLuminance ; 
    % producing BW for bars in either position.  0 inside the patch, 255 outside.
    BW(:,:,2) = 255-255*PMB.polymask{c,f,p} ;
    if PMB.edgeSmoothing
        BWtemp = padarray(BW(:,:,2),[100 100], 255, 'both');
        BWtemp = imfilter(BWtemp, h);
        BW(:,:,2) = BWtemp(101:end-100, 101:end-100);
    end

    PMB.masktex(c,f) = Screen('MakeTexture', screen, BW ) ;

%     if PMB.save_screen==1
%         ScreenImg = DrawScreen(screenRect,[8,6]); %screenRect=[0,0,1600,900]
%         ScreenImg(:,:,1) = ScreenImg(:,:,1) - PMB.polymask{c,f,p} ;
%         ScreenImg(:,:,3) = ScreenImg(:,:,1);
%         imwrite(ScreenImg,[save_dir_maskscreen 'screenmask_config' num2str(c) '.png'],'png');
%     end
    if Param.Dichoptic
            ScreenImg_merge(:,:,2,f) = ScreenImg_merge(:,:,2,f) - PMB.polymask{c,f,p} ;
            ScreenImg_merge(:,:,3,f) = ScreenImg_merge(:,:,2,f);
    else
        ScreenImg_merge(:,:,2,1) = ScreenImg_merge(:,:,2,1) - PMB.polymask{c,f,p} ;
        ScreenImg_merge(:,:,3,1) = ScreenImg_merge(:,:,2,1);
    end
%     end

end
end

end % for p = 1:PMB.positions(1)

if PMB.save_screen==1
    for i=1:length(screen)
        imwrite(ScreenImg_merge(:,:,:,i),[save_dir_maskscreen 'screenmask' num2str(i) '_merge.png'],'png');
    end
end
% if show_screen_only==1
%     imshow(ScreenImg_merge);
%     return
% end

PMB.seqbars=[]; seqbarsdichtmp=[]; seqbarsdich=[]; rowstokeep=[];
[A,B]=meshgrid( 1:PMB.locations , 1:PMB.locations );
seqbarstmp = reshape( cat(2,A,B) ,[],2 ); % sequences of dichoptic phi motion bars and "pseudostatic" bars
for ii=1:size(seqbarstmp,1)
    if seqbarstmp(ii,1)~=seqbarstmp(ii,2); rowstokeep = [rowstokeep,ii]; end;
end; seqbarstmp = seqbarstmp(rowstokeep,:);
% if Param.Dichoptic
% %     additionalStaticBars = repmat([1,1; 2,2; 3,3] ,3,1);
    PseudoStaticBars = 1:PMB.locations; PseudoStaticBars = [PseudoStaticBars;PseudoStaticBars]';
%     seqbarstmp=cat(1,seqbarstmp,additionalStaticBars);
% end
seqbarsdichtmp = [ PseudoStaticBars; seqbarstmp];
seqbarsdichtmp = [seqbarsdichtmp ; seqbarsdichtmp]; % duplicate for screensequence [1 2] and [2 1]
for ii = 1:numel(PMB.rotation_deg)
    seqbarsdich = [seqbarsdich ; seqbarsdichtmp+(ii-1)*PMB.locations ] ;
end
seqbarsmono = seqbarsdich;
PMB.seqbars = [seqbarsmono; seqbarsdich];
if     PMB.Black == 1 && PMB.White == 1
    PMB.seqbars = [PMB.seqbars ; PMB.seqbars ] ; % duplicate for black and white.
end
if     PMB.Black == 1 && PMB.White == 1
    PMB.seqcolor( 1:PMB.nbr_stimuli/2 ,1) = PMB.BlackLuminance ;       % first half black
    PMB.seqcolor( PMB.nbr_stimuli/2 +1 : PMB.nbr_stimuli ,1) = PMB.WhiteLuminance ;   % second half white
elseif PMB.Black == 1 && PMB.White == 0
    PMB.seqcolor( 1:PMB.nbr_stimuli ,1) = PMB.BlackLuminance ;
elseif PMB.Black == 0 && PMB.White == 1
    PMB.seqcolor( 1:PMB.nbr_stimuli ,1) = PMB.WhiteLuminance ;
end
if Param.Dichoptic
    screensequencedich = [ repmat([1 2], length(seqbarsdichtmp)/2,1); repmat([2 1], length(seqbarsdichtmp)/2,1)];
    screensequencemono = [screensequencedich(:,1),screensequencedich(:,1)]; screensequencemono([1:PMB.locations, PMB.locations^2+1:PMB.locations^2+PMB.locations],2)=0;
    screensequencedich = repmat( screensequencedich, numel(PMB.rotation_deg),1);
    screensequencemono = repmat( screensequencemono, numel(PMB.rotation_deg),1);
%     screensequencemono = screensequencedich; screensequencemono(:,2)=0;
    PMB.screensequence = [ screensequencedich; screensequencemono];
%     PMB.screensequence = repmat(PMB.screensequence, (PMB.Black+PMB.White)*numel(PMB.rotation_deg),1);
else
    PMB.screensequence = [1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];
    PMB.screensequence = repmat(PMB.screensequence, (PMB.Black+PMB.White)*numel(PMB.rotation_deg),1);
end
for r = 1:PMB.n_reps
    % seqpositionstmp is [1 1; 1 2; 1 3; 1 4; 2 1; ...; 6 4]:
    switch Param.seqmode
        case 'sequential'
            PMB.stimseq(r,:) = 1 : PMB.nbr_stimuli ;
        case 'random'
            PMB.stimseq(r,:) = randperm( PMB.nbr_stimuli ) ;
    end     
end

Param.stimSeq(end+1,1) = { repmat(PMB.stim_id, 1,PMB.nbr_stimuli*PMB.n_reps) };
% stimSeqtmp4 = [ repmat(PMB.stim_id, 1,PMB.nbr_stimuli*PMB.n_reps) ];
PMB.cnt=1;  

if PMB.save_screen==1
    if show_screen_only==1
        Screen('CloseAll');
        Priority(0);
        ShowCursor;
        hwb = waitbar(0,'Saving screen configurations...');
    end
    for s=1:PMB.nbr_stimuli/(PMB.Black+PMB.White)
        if show_screen_only==1
            waitbar(s/(PMB.nbr_stimuli/(PMB.Black+PMB.White)), hwb);
        end
        ScreenImg = DrawScreen(Param.screenRect,patches); %screenRect=[0,0,1600,900]
        for f = [1,2]%PMB.screensequence % [1,2]
            if     f==1%PMB.screensequence(1) % first bar in light green
            ScreenImg(:,:,1) = ScreenImg(:,:,1) - PMB.polymask{PMB.seqbars(s,f),f,1} ;
            elseif f==2%PMB.screensequence(2) % second bar in dark green
            ScreenImg(:,:,1) = ScreenImg(:,:,1) - 0.7*PMB.polymask{PMB.seqbars(s,f),f,1} ;
            ScreenImg(:,:,2) = ScreenImg(:,:,2) - 0.4*PMB.polymask{PMB.seqbars(s,f),f,1} ;
            end
            ScreenImg(:,:,3) = ScreenImg(:,:,1);
            imwrite(ScreenImg,[save_dir_maskscreen 'screenmask_config' num2str(s) '.png'],'png');
        end
    end
    if show_screen_only==1; close(hwb); end;
end

if exist('scriptNamePhi','var')
    BackupStimulusScript( scriptPathPhi, scriptNamePhi, save_dir )
end
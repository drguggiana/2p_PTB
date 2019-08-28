function [ RFM, Param ] = SetRFMapping( screen, RFM, Param )

    [scriptPath, scriptName] = fileparts(mfilename('fullpath'));


    RFM.frame_patch = round(RFM.patch_time/Param.ifi); % # frames per patch.
    RFM.frame_interpatch = round(RFM.interpatch_time/Param.ifi);
    % to position the stimulation in the screen. Center of the screen:
    RFM.mouseCenter(1,:) = [ (Param.screenRect(1,3)-Param.screenRect(1,1))/2  (Param.screenRect(1,4)-Param.screenRect(1,2))/2];%19/25 ] ; % in pixel coordinates (position the mouse pointer on the screen an use GetMouse in MatLab)
    if     strcmp(inputname(2),'RFM2')
        RFM.mouse_dist_px(1) = round(Param.mouse_dist_cm(2)/Param.pixel_cm);
    elseif strcmp(inputname(2),'RFM')
        RFM.mouse_dist_px(1) = round(Param.mouse_dist_cm(1)/Param.pixel_cm);
    end

    if RFM.save_screen==1
        save_dir_maskscreen = [Param.save_dir filesep 'screen_masks_' Param.datetime_suffix filesep];
        if ~isdir(save_dir_maskscreen), mkdir(save_dir_maskscreen); end;
    end

    
    
    
    patches = [18 6];
    for i=1:length(screen)
        ScreenImg_merge(:,:,:,i) = DrawScreen(Param.screenRect, patches);
    end

    % calculate the patch shapes

    [lo,la]=patches_deg(RFM.n_patches, RFM.fov, RFM.view_offset - RFM.fov/2 , RFM.rel_patch_size);  %Ale: produces vectors with all corner points for the patches

    switch RFM.gnomonic
        case 1
            [x,y]=  pr_gnomonic(reshape(lo, [],1),reshape(la, [],1));   %Ale: pr_gnomonic: coordinates with the Gnomonic non conformal projection;                                               %     reshape lo and la (that are matrices #tot_patches-by-4) in one single column
            xx = reshape(x,[],4);   %Ale: restore the matrix with 4 columns
            yy = reshape(y,[],4);
        case 0
            xx(:,:) = lo;
            yy(:,:) = la;
    end

    RFM.masktex=zeros(1,prod(RFM.n_patches));
    % BW is a matrix that contains the white and black patches (BW(:,:,1)==white, BW(:,:,2) == black)
    % save BW and load every time will speed up
    %     BW1=zeros(screenRect(4),screenRect(3)-screenRect(1),2,prod(n_patches));

    for p = 1:prod(RFM.n_patches)

        RFM.polymask{p,1} = poly2mask(xx(p,:).* RFM.mouse_dist_px(1) + RFM.mouseCenter(1,1), yy(p,:).* RFM.mouse_dist_px(1) + RFM.mouseCenter(1,2), Param.screenRect(1,4),Param.screenRect(1,3)-Param.screenRect(1,1));
        % BW are LA planes where A is alpha, the transparency of a pixel.
        % Alpha values range between zero(=fully transparent) and 255(=fully
        % opaque).
        % polymask: 1 inside the bar, 0 outside.
        % Multiplying by mouse_dist_px you take into account the distance from
        % the screen and rescale all the measures in pixels to obtain the
        % proper values in degrees.
        BW = zeros(Param.screenRect(1,4),Param.screenRect(1,3)-Param.screenRect(1,1),2);
        BW(:,:,2) = 255-255*RFM.polymask{p,1};
        %BW(:,:,2) = 0 inside the patch, 255 outside.
        BW(:,:,1) = BW(:,:,1) + RFM.BackgroundLuminance;
    % %         BW2(:,:,p)=BW(:,:,2);
    %         BW1(:,:,p)=BW(:,:,1);

        RFM.masktex(p) = Screen('MakeTexture', screen(1), BW);

            ScreenImg_merge(:,:,2,1) = ScreenImg_merge(:,:,2,1) - RFM.polymask{p,1} ;
            ScreenImg_merge(:,:,3,1) = ScreenImg_merge(:,:,2,1);

    end

    RFM.ScreenImg_merge = ScreenImg_merge;
    if RFM.save_screen==1
        for i=1:length(screen)
            imwrite(ScreenImg_merge(:,:,:,i),[save_dir_maskscreen 'screenmask' num2str(i) '_merge.png'],'png');
        end
    end

    j=1;
    while j<=RFM.n_reps
        stimseq(j,:)= randperm(prod(RFM.n_patches)*2);  %Ale: random sequence of prod(n_patches)*2 integers.
    %     stimseq(j,:)= 1:prod(n_patches)*2 ;
        % check if consecutive frames are    %Ale: what are you checking??
        check_aux=stimseq(j,:);
        aux=(stimseq(j,:)<=prod(RFM.n_patches));   %Ale: vector of 0 and 1.
        check_aux(aux)=stimseq(j,aux)+prod(RFM.n_patches);   %Ale: random sequence of prod(n_patches)*2 integers between prod(n_patches)+1 and prod(n_patches)*2.
        if (sum(check_aux(1:2:end)==check_aux(2:2:end))>0)|( sum(check_aux(2:2:end-1)==check_aux(3:2:end-1))>0);    %Ale: why this?
            stimseq(j,:)=zeros(1,(prod(RFM.n_patches)*2));
            j;
        else j=j+1;
            RFM.stimseq(j-1,:) = stimseq(j-1,:);
        end

    end
    % the sequence of patches is columnwise! %%%%
    RFM.patchseq=[1:(prod(RFM.n_patches)),1:(prod(RFM.n_patches))];  % vector [1 2 3... prod(n_patches) 1 2 3... prod(n_patches)].
    RFM.colorseq(1:prod(RFM.n_patches)) = RFM.BlackLuminance;                          % first prod(n_patches) black
    RFM.colorseq(prod(RFM.n_patches)+1:prod(RFM.n_patches)*2) = RFM.WhiteLuminance;   % second prod(n_patches) white

    Param.stimSeq(end+1,1) = { repmat(RFM.stim_id, 1,prod(RFM.n_patches)*2*RFM.n_reps) };
    RFM.cnt=1; 


    if exist('scriptName','var')
        if (length(Param.StimProtocols.RFMapping)>1 && strcmp(inputname(1),'RFM'))...
            || length(Param.StimProtocols.RFMapping)==1
            save_dir = evalin('caller', 'save_dir');
            BackupStimulusScript( scriptPath, scriptName, save_dir )
        end
    end

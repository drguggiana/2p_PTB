    
function DG = SetMask(eye, DG, Param, pixperdeg, win)
% Draw smoothed circular mask.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Only in case of DGD, be careful to set:
%   eye=1 -> left screen
%   eye=2 -> right screen.
% For DG and DG2, just use eye=1.

fsize = 31;

DG.maskSize_px(eye)   = DG.maskSize_deg(eye) * pixperdeg ;
DG.maskRadius_px(eye) = round(DG.maskSize_px(eye)/2);
DG.maskSize_px(eye)   = 2*DG.maskRadius_px(eye) +1 ;

tmasksize = 2*DG.maskSize_px(eye);
% calculate offset of grating in pixels:
DG.maskCenter_offset_px(eye,:) = DG.maskCenter_offset(eye,:) .* [Param.screenRect(3), Param.screenRect(4)] ;

% Important: gratingsize and gratingRect have to be the same size,
% otherwise the grating will be rescaled
DG.gratingsize(eye,:) = [tmasksize, tmasksize];


% Smoothed circular aperture mask

% Create transparency mask (enlarged because we will smooth the
% edge of the actual circular aperture):
mask = ones(tmasksize, tmasksize, 2) * DG.BackgroundLuminance;
% Circle with radius DG.maskRadius_px centered at
% (tmasksize/2,tmasksize/2) in image tmasksizeXtmasksize:
[rr, cc] = meshgrid(1:tmasksize);
C = sqrt((rr-tmasksize/2).^2 + (cc-tmasksize/2).^2) <= DG.maskRadius_px(eye) ;
Csm = imfilter(double(C), fspecial('average',[fsize fsize]));
% Final transparency mask (0=transparent, 255=oapque):
mask(:,:,2) = 255 * (1 - Csm);
%         figure; imagesc(mask(:,:,2));

DG.masktex{eye} = Screen('MakeTexture', win, mask);

% Definition of the drawn rectangle on the screen, centered on the
% screen center plus maskCenter_offset:
[DG.gratingRect(eye,:),dh,dv] = CenterRect([0 0 tmasksize tmasksize], Param.screenRect);
DG.gratingRect(eye,:) = DG.gratingRect(eye,:) + [DG.maskCenter_offset_px(eye,1) DG.maskCenter_offset_px(eye,2) DG.maskCenter_offset_px(eye,1) DG.maskCenter_offset_px(eye,2)];
	

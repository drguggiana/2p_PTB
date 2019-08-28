function [X,Y]=pr_gnomonic(Long,Lat,R,CenterVec);
%--------------------------------------------------------
% pr_gnomonic function    Project coordinates with the
%                     Gnomonic non conformal projection
%                     for set of longitude and latitude.
%                     This is a nonconformal projection
%                     from a sphere center in which
%                     orthodromes are stright lines.
% Input  : - vector of Longitude, in radians.
%          - vector of Latitude, in radians.
%          - Scale radius, default is 1.
%          - central coordinate vector [Long_center,Lat_center]
% Output : - vector of X position
%          - vector of Y position 
%    By : Eran O. Ofek        July 1999
%--------------------------------------------------------
if (nargin==4),
   % no default
elseif (nargin==3),
   CenterVec = [0 0];
elseif (nargin==2),
   CenterVec = [0 0];
   R = 1;
else
   error('Illigal number of argument');
end
%CenterVec = [-0.2514 0.4727];
% extract central coordinates.

% Change CenterVec by data middle
%--------------------------------
%CenterVec(1) = (max(Long)+min(Long)).*0.5;
%CenterVec(2) = (max(Lat)+min(Lat)).*0.5;
%--------------------------------

Long1 = CenterVec(1);
Lat1  = CenterVec(2);
% R is really R.*S (R-radius, S-scale factor)
CosC = sin(Lat1).*sin(Lat) + cos(Lat1).*cos(Lat).*cos(Long-Long1);
%X = cos(Lat).*sin(Long-Long1)./CosC;
X = sin(Long-Long1)./CosC;
% Y = sin(Long)./CosC;
% Y = sin(Lat).*sin(Long-Long1)./CosC;
%X = sin(Lat)./CosC;
Y = (cos(Lat1).*sin(Lat) - sin(Lat1).*cos(Lat).*cos(Long-Long1))./CosC;
%X = sin(Lat).*cos(Long)./cos(Lat);
%Y = sin(Lat).*sin(Long)./cos(Lat);



function ScreenImg = DrawScreen(screenRect, positions)

if nargin==1
    positions=[6 4];
elseif nargin==0
    screenRect=[0,0,1600,900];
    positions=[6 4];
end

size_grating = [ fix(screenRect(3)/positions(1))  fix(screenRect(4)/positions(2)) ] ;

for py = 1 : positions(2)
for px = 1 : positions(1)
	pos_coord(px,py,:) = [ screenRect(1) + (px-1)*size_grating(1) ...
						   screenRect(2) + (py-1)*size_grating(2) ...
					       screenRect(1) + (px)  *size_grating(1) ...
					       screenRect(2) + (py)  *size_grating(2) ] ;
end
end
ScreenImg = ones(screenRect(4) , screenRect(3) , 3); %screenRect=[0,0,1600,900]
% ScreenImg( pos_coord(1,2:end,2) , : ,:) = 0; %y
% ScreenImg( : , pos_coord(2:end,1,1) ,:) = 0; %x
width=5; % line width will be width*2+1
vert=pos_coord(2:end,1,1);
horiz=pos_coord(1,2:end,2);
for v=1:length(vert)
    ScreenImg( :, vert(v)-width:vert(v)+width, :) = 0;
end
for h=1:length(horiz)
    ScreenImg( horiz(h)-width:horiz(h)+width, :, :) = 0;
end
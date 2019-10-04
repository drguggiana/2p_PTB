function masktex = CreateCircularAperature(radius, win, Param)

    white = 255;
    grey = 127;
    screenYpix = Param.screenRect(4) / 2;
    
    % Make the mask to cover the texture
    rad = radius * 2 * pi;
    [xm, ym] = meshgrid(-(screenYpix-1):screenYpix, -(screenYpix-1):screenYpix);
    [s1, s2] = size(xm);
    
    % Define the circular aperature as a boolean mask
    circle = xm.^2 + ym.^2 <= rad^2;
    
    mask = ones(s1, s2, 2)*grey;
    mask(:,:,2) = ~circle .* white; 
    masktex = Screen('MakeTexture', win, mask);        
end
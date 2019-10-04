function masktex = CreateCircularAperature(radius, win, xm, ym)

    white = 255;
    grey = 127;
    
    % Make the mask to cover the texture
    rad = radius * 2 * pi;
    [s1, s2] = size(xm);
    
    % Define the circular aperature as a boolean mask
    circle = (xm.^2 + ym.^2) <= rad^2;
    
    mask = ones(s1, s2, 2)*grey;
    mask(:,:,2) = ~circle .* white; 
    masktex = Screen('MakeTexture', win, mask);        
end
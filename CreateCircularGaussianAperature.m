function [masktex, fullWindowMask] = CreateCircularGaussianAperature(radius, win, Param)

    white = 255;
    grey = 127;
    screenYpix = Param.screenRect(4) / 2;
    
    % Make the mask to cover the texture
    rad = 2*pi*radius;
    [xm, ym] = meshgrid(-screenYpix:screenYpix, -screenYpix:screenYpix);
    [s1, s2] = size(xm);
    
    % Define the circular aperature as a boolean mask
    circle = xm.^2 + ym.^2 <= rad^2;
    
    mask = ones(s1, s2, 2)*grey;
    mask(:,:,2) = ~circle .* white; 
    masktex = Screen('MakeTexture', win, mask);

    % Make the full window texture
    fullWindowMask = Screen('MakeTexture', win, ...
        ones(Param.screenRect(3), Param.screenRect(4)) .* grey);
        
end
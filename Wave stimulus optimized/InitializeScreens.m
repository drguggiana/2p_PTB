function [w, screenNumber, colors, ifi, oldRes] = InitializeScreens( screenNumber, PrefRes, SphereCorrection )
    
    % Check ig OpenGL can be used
    AssertOpenGL;

    % Find which screen to project on
    screens=Screen('Screens');
    if strcmpi(screenNumber,'Max')
        screenNumber = max(screens);
    end
    
%     % Set resolution
%     oldRes = SetResolution(screenNumber,PrefRes(1),PrefRes(2));
    oldRes = [PrefRes(1),PrefRes(2)];
     
    % Get color range
    colors.white=WhiteIndex(screenNumber);
    colors.black=BlackIndex(screenNumber);

    % Open on-screen window
    if isempty(SphereCorrection)
        
        % Open flat screen
        w = Screen('OpenWindow',screenNumber, 0,[],32,2);
        
    else
        
        % Open Sphere corrected screen
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'AllViews', 'GeometryCorrection', SphereCorrection);
%         w = PsychImaging('OpenWindow',screenNumber, 0,[],32,2);
        w = PsychImaging('OpenWindow',screenNumber, 0);
        
    end
    
%     Screen('LoadNormalizedGammaTable',w,linspace(0,1,256)'*ones(1,3));
%     Clut = [(254:-1:-1)*0; -1:254; -1:254]' * 257;
%     BitsPlusSetClut(w,Clut);    
    
    ifi = Screen('GetFlipInterval', w );

    
    Screen('FillRect',w, colors.black);
    Screen('Flip', w);
    Screen('FillRect',w, colors.black);
    Screen('TextFont', w, 'Arial');
    Screen('TextSize', w, 20);
end

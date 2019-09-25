% This function differs from the normal PresentLooming.m by drawing an
% expanding transparent circle over a checkered background to create a
% an expanding looming checkerboard.


function [stimtype, vbl] = PresentLoomingRadialCheckerboard( nidaq, screen, vbl, Param, stimtype, tar_eye )

    [scriptPath, scriptName] = fileparts(mfilename('fullpath'));
    tic;

    if     Param.Dichoptic && Param.stereoMode==0 
    %     if     strcmp(inputname(5),'DG') || strcmp(inputname(5),'DG3')
            win=screen(tar_eye);
    %     elseif strcmp(inputname(5),'DG2')|| strcmp(inputname(5),'DG4')
    %         win=screen(2);
    %     end
        winStereo = 0; % for Photodiode
        Screen('FillRect', screen(1), stimtype.BackgroundLuminance);
        Screen('FillRect', screen(2), stimtype.BackgroundLuminance);
        vbl = Screen('Flip', screen(1));
        vbl = Screen('Flip', screen(2));
    % elseif Param.Dichoptic && ismember(Param.stereoMode,[1 4])
    %     if     strcmp(inputname(5),'DG')  || ...
    %            strcmp(inputname(5),'DG3') || ...
    %            strcmp(inputname(5),'DG5')
    %         winStereo = 0;
    %     elseif strcmp(inputname(5),'DG2') || ...
    %            strcmp(inputname(5),'DG4') || ...
    %            strcmp(inputname(5),'DG6')
    %         winStereo = 1;
    %     end
    %     win=screen(1);
    %     Screen('SelectStereoDrawBuffer', win, 0);
    %     Screen('FillRect', win, stimtype.BackgroundLuminance);
    %     Screen('SelectStereoDrawBuffer', win, 1);
    %     Screen('FillRect', win, stimtype.BackgroundLuminance);
    %     vbl = Screen('Flip', win);
    else
        win=screen(1);
        winStereo = 0;
        vbl = Screen('Flip', win);
    end

    % if stimtype.drawMask == 0
    % 	Screen('BlendFunction', win, GL_ONE, GL_ZERO);
    % else % if gaussian transparency mask:
    % 	% Enable alpha blending for combination of the gaussian aperture
    % 	% with the drifting grating:
    % 	Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    % end

    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    % Screen('BlendFunction', win, GL_ONE, GL_ZERO);

    % Disable alpha-blending, restrict following drawing to alpha channel:
    % Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]);

    
    % draw the background checkerboard
    % Draw our texture to the screen
    Screen('DrawTexture', win, LO.radialChecker);
    
    %define the rectangle that will bind the circle, using its max size
    baseRect = [0 0 1 1];
    [xCenter, yCenter] = RectCenter(stimtype.gratingRect);

    % Center the rectangle on the centre of the screen
    centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);

    %define the max Diameter
    %max_diam = max(screenYpixels);

    % Sync us and get a time stamp
    vbl = Screen('Flip', win);

    %define the coordinates of the circle center
    circle_center_x = 0.5*stimtype.gratingsize(1);
    circle_center_y = 0.5*stimtype.gratingsize(2);

    % %define the final size of the circle
    % max_size = 0.33*screenXpixels;

    % %define the scale factor in pixels/frame
    % scaleFactor = max_size/n_frames;

    test=Param.test;
    tempo=[0 0];

    %transpose the times and speeds arrays
    seqtimes = stimtype.seqtimes';
    seqspeeds = stimtype.seqspeeds';
    seqspeedorder = stimtype.paramorder(:,:,1)';

    vblendtime = vbl + seqtimes(stimtype.cnt) - Param.ifi;

    %define the scale factor for the expanding circle
    scaleFactor = seqspeeds(stimtype.cnt);

    %get the current color - this is fully transparent
    stim_color = [0,0,0,0];

    %also print the current expansion speed
    fprintf(strcat('Current exp speed:',...
        num2str(stimtype.expansion_speeds(seqspeedorder(stimtype.cnt))),...
        '/Current color:',num2str(stimtype.seqcolors(stimtype.cnt,:)),...
        '\r\n'))
    %         for i = 1 : frame_stimulus_dur
    frameNr = 0;

    while vbl < vblendtime
        frameNr = frameNr + 1;
       
        % variable intensity circle
        Screen('FillOval', win, stim_color, centeredRect);

        %update the size of the circle
        baseRect(3:4) = baseRect(3:4) + scaleFactor;

        % Center the rectangle on the centre of the screen
        centeredRect = CenterRectOnPointd(baseRect, circle_center_x, circle_center_y);

        if Param.usePhotodiode
            if ismember(Param.stereoMode,[1 4])
                Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
            end
            Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
        end

        Screen('DrawingFinished', win);

        if frameNr == 1
            StimSignal( Param, test, nidaq, stimtype.stim_id )
        end

        % Show it at next retrace:
        vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);

        Exit_If_Shift(test,nidaq);
    end
    tempo(1) = toc;

    StimSignal( Param, test, nidaq, 0 )

    % Post-stim:
    tic
    Screen('FillRect', win, stimtype.BackgroundLuminance);
    while  tempo(2) < stimtype.poststim_time
        tempo(2) = toc;
        Exit_If_Shift(test,nidaq);
        vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);  
    end
    stimtype.tempo{stimtype.cnt}=tempo;

    if exist('scriptName','var') && stimtype.cnt==1
        save_dir = evalin('caller', 'save_dir');
        BackupStimulusScript( scriptPath, scriptName, save_dir )
    end
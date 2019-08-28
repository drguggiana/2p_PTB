function RunWave(Settings,StimSettings)

%     OutputDir = [Settings.ExptDir filesep Settings.ExptName];
%     x = dir(OutputDir);
%     if length(x) == 0
%         mkdir(OutputDir);
%     else
%         MatFiles = dir([OutputDir filesep '*.mat']);
%         if ~isempty(MatFiles)
%             disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
%             disp('Detected possible data-files in experiment directory:');
%             for f = 1:length(MatFiles)
%                 disp([' * ' OutputDir filesep MatFiles(f).name]);
%             end
%             Answer = input('Can these be overwritten (1=yes/0=no): ');
%             if Answer ~= 1
%                 disp('Aborted.. please restart for a new experiment');
%                 return
%             end
%         end
%     end

    % Stimulus settings
    Type = StimSettings.WaveType;
    BarAperture = StimSettings.BarSize;
    BarSpeed = StimSettings.BarSpeed;
    StimulusRepeat = Settings.NumTrials;

    Options.RepWait = StimSettings.LoopDelay;
    Options.shape = StimSettings.GratShape;
    Options.spatialF = StimSettings.SpatFreq;
    Options.temporalF = StimSettings.TempFreq;
%     Angles = StimSettings.Angles;
    
    NumFlickerFrames = 60/StimSettings.FlickerFreq;
%     NumDirectionFrames = 60/StimSettings.DirAltFreq;
    
    Options.contrastmin = StimSettings.MinInt;
    Options.contrastmax = StimSettings.MaxInt;
    BackgroundColor = StimSettings.BGint;
    
    Options.scrx = StimSettings.ScrResX;
    Options.scry = StimSettings.ScrResY;
    Options.scrwidth = StimSettings.ScrSizeX;
    Options.scrheight = StimSettings.ScrSizeY;
    Options.scrdist = StimSettings.ScreenDist;
    MarginX = round((StimSettings.ScrResX-StimSettings.PresAreaX)/2);
    MarginY = round((StimSettings.ScrResY-StimSettings.PresAreaY)/2);
    Options.xsize = round(Options.scrx*1.5);
    Options.ysize = round(Options.scry*1.5);
    
    Options.ArdPinNr = Settings.ArdPinNr;

    save(Settings.stim_fname, 'Settings', 'StimSettings', 'Options' );


    if Settings.DipslayBlueGreen == 1
        BackgroundColor = [0 1 1] * BackgroundColor;
    end

    % Screen & pixel size in degrees
    ScrWidthDegree  = 2 * atan( (0.5*Options.scrwidth) /Options.scrdist) * (180/pi);
    ScrHeightDegree = 2 * atan( (0.5*Options.scrheight)/Options.scrdist) * (180/pi);
    PixelWidthDegree  = ScrWidthDegree  / Options.scrx;
    PixelHeightDegree = ScrHeightDegree / Options.scry;

    % Bar speed in pixels
    HorzBarSpeedPix = (BarSpeed / PixelHeightDegree) / 60; % 60Hz screen
    VertBarSpeedPix = (BarSpeed / PixelWidthDegree)  / 60;

    % Bar height in pixels
    HorzBarHeightPix = round( BarAperture / PixelHeightDegree );
    VertBarWidthPix  = round( BarAperture / PixelWidthDegree  );

    % Steps in pixels for bar movement
    HorzBarYsteps = round( (-HorzBarHeightPix+MarginY):HorzBarSpeedPix:(Options.scry-MarginY) );
    VertBarXsteps = round( (-VertBarWidthPix-MarginX):VertBarSpeedPix:(Options.scrx+MarginX));
    NumHorzBarFrames = length(HorzBarYsteps);
    NumVertBarFrames = length(VertBarXsteps);

    % Set indices for horizontal and vertical bars
    NumMovieFramesHor = (StimulusRepeat*NumHorzBarFrames);
    MovieDurationHor = NumMovieFramesHor/60;
    
    NumMovieFramesVer = (StimulusRepeat*NumVertBarFrames);
    MovieDurationVer = NumMovieFramesVer/60;

    NumMovieFramesBoth = (StimulusRepeat*max([NumHorzBarFrames NumVertBarFrames]));
    MovieDurationBoth = NumMovieFramesBoth/60;
    
    HorFrameIx = mod( (1:NumMovieFramesBoth)-1, NumHorzBarFrames ) + 1;
    VerFrameIx = mod( (1:NumMovieFramesBoth)-1, NumVertBarFrames ) + 1;
    HorFrameIxRev = abs(HorFrameIx-NumHorzBarFrames)+1;
    VerFrameIxRev = abs(VerFrameIx-NumVertBarFrames)+1;
    
    
    if strcmpi( Type, 'Bar' ) == 1
    end
    
    if strcmpi( Type, 'Flicker' ) == 1 || strcmpi( Type, 'FlickerShifty' ) == 1
        % get horizontal and vertical grating
        HorGrat = SimpleGrating( 180, Options );
        VerGrat = SimpleGrating( 90, Options );
        Flicker1 = HorGrat + VerGrat;
        Flicker1(Flicker1~=Options.contrastmax) = Options.contrastmin;
        Flicker2 = abs(Flicker1-Options.contrastmax);
        if Settings.DipslayBlueGreen == 1
            Flicker1 = repmat(Flicker1,[1 1 3]);
            Flicker1(:,:,1) = 0;
            Flicker2 = repmat(Flicker2,[1 1 3]);
            Flicker2(:,:,1) = 0;
        end
        FlickerFrameIndices = mod( floor(((1:NumMovieFramesBoth)-1) ./ NumFlickerFrames), 2 ) + 1;
        RandShiftFrameIndices = ceil( (1:NumMovieFramesBoth) ./ (2*NumFlickerFrames));
    end
    
    if strcmpi( Type, 'MovGrat' ) == 1
        % loop for angle of the grating
        for a = 1:length(Angles)
            [ Grating{a}, xSteps{a}, ySteps{a}, Specs ] = SimpleGrating( Angles(a), Options );
            if Settings.DipslayBlueGreen == 1
                Grating{a} = repmat(Grating{a},[1 1 3]);
                Grating{a}(:,:,1) = 0;
            end
        end
        MovieFrameIndices = mod( (1:NumMovieFramesBoth)-1, length(xSteps{1}) ) + 1;
        AngleFrameIndices = mod( ceil( (1:NumMovieFramesBoth)/NumDirectionFrames )-1, length(Angles) )+1;
    end
    
    
    %% start psychophysics toolbox
    try
        
        AssertOpenGL;
        Screen('Preference', 'SkipSyncTests', 1);
        [w, screenNumber, colors, ifi, oldRes] = InitializeScreens( 'Max', ...
            [Options.scrx Options.scry], Settings.CurvCorrectFile );
        AssertGLSL;
        if ~isempty(Settings.GammaCorrectionTable)
            Screen('LoadNormalizedGammaTable', w, Settings.GammaCorrectionTable*[1 1 1]);
        else
            disp(' ');
            disp('WARNING: No gamma correction implemented!!');
            disp(' ');
        end
%         disp('Note: Acquisition set for monochrome camera, e.g. DALSA 1M60');
        Screen('FillRect',w, 0);
        Screen('Flip', w); 

        %% load the stimulus into video memory (tex)
        if strcmpi( Type, 'Bar' ) == 1
            if Settings.DipslayBlueGreen == 1
                ScrFill = ones(Options.scry,Options.scrx,3)*Options.contrastmax;
                ScrFill(:,:,1) = 0;
                Tex = Screen('MakeTexture', w, ScrFill);
            else
                Tex = Screen('MakeTexture', w, ones(Options.scry,Options.scrx)*Options.contrastmax);
            end
        elseif strcmpi( Type, 'Flicker' ) == 1 || strcmpi( Type, 'FlickerShifty' ) == 1
            Tex(1) = Screen('MakeTexture', w, Flicker1);
            Tex(2) = Screen('MakeTexture', w, Flicker2);
        elseif strcmpi( Type, 'MovGrat' ) == 1
            for a = 1:length(Angles)
                Tex(a) = Screen('MakeTexture', w, Grating{a});
            end
        end
        
        % set the screen to a single background color
        Screen('FillRect',w, BackgroundColor);
        Screen('Flip', w);
        
        disp('Press any key to start animation')
        while ~KbCheck
        end
        
        disp(['Wait ' num2str(StimSettings.prestimInterval,'%2.0f') ' seconds to recover from light flash of psychtoolbox']);
        tic;
        while toc < StimSettings.prestimInterval
            Screen('FillRect',w, BackgroundColor);
            Screen('Flip', w);
            pause(0.1);
        end
        
        %% Start acquisition loop
                
        % Run stimuli and acquire
        if StimSettings.XYsimult == 1
            RunIDs = [5,6];
            RunSpecs = {'Horizontal down & Vertical right','Horizontal up & Vertical left'};
            RunName = {'HorDownVertRight','HorUpVertLeft'};
        else
%             disp('WARNING: RUNNING 1 DIRECTION ONLY >> TESTING!!!');
%             RunIDs = [1];
%             RunSpecs = {'Horizontal down'};
%             RunName = {'HorDown'};
            RunIDs = [1,3,2,4];
            RunSpecs = {'Horizontal down','Vertical right','Horizontal up','Vertical left'};
            RunName = {'HorDown','VertRight','HorUp','VertLeft'};
        end
        
        
        % Loop stimulus runs
        for r = 1:length(RunIDs)
           
            % Calibrate imaging device to acquire X frames and average them
            if RunIDs(r) == 1 || RunIDs(r) == 2
                NumMovieFrames = NumMovieFramesHor;
            elseif RunIDs(r) == 3 || RunIDs(r) == 4
                NumMovieFrames = NumMovieFramesVer;
            elseif RunIDs(r) == 5 || RunIDs(r) == 6
                NumMovieFrames = NumMovieFramesBoth;
            end
            
            if r > 1
                disp('Wait 10 seconds to recover from light flashes');
                tic;
                while toc < 10
                    Screen('FillRect',w, BackgroundColor);
                    Screen('Flip', w);
                    pause(0.1);
                end
            end
            
%             % Image acquisition settings
%             set(Settings.CamId, 'ROIPosition', Settings.ROI );
%             set(Settings.CamId, 'FramesPerTrigger', 1);
%             set(Settings.CamId, 'TriggerRepeat',inf);
%             Settings.CamId.ReturnedColorspace = 'grayscale';
%             imaqmem(25000000000);

            % Now start the actual imaging
            disp(['- Start imaging. Stimulus: ' RunSpecs{r}]);
%             start(Settings.CamId); 
%             tic;
            StimSignal( Settings, Settings.test, [], 1 )
            
            % Show the stimulus
            if strcmpi( Type, 'Bar' ) == 1
                [LoopTimestampHor,LoopTimestampVer] = ...
                    BarBar( RunIDs(r), w, ifi, Tex, NumMovieFrames, BackgroundColor, ...
                        HorzBarHeightPix, VertBarWidthPix, HorzBarYsteps, VertBarXsteps, ...
                        HorFrameIx, HorFrameIxRev, VerFrameIx, VerFrameIxRev, Options );
            elseif strcmpi( Type, 'Flicker' ) == 1
                [LoopTimestampHor,LoopTimestampVer] = ...
                    FlickerBar( RunIDs(r), w, ifi, Tex, NumMovieFrames, FlickerFrameIndices, BackgroundColor, ...
                        HorzBarHeightPix, VertBarWidthPix, HorzBarYsteps, VertBarXsteps, ...
                        HorFrameIx, HorFrameIxRev, VerFrameIx, VerFrameIxRev, Options );
            elseif strcmpi( Type, 'FlickerShifty' ) == 1
                [LoopTimestampHor,LoopTimestampVer] = ...
                    FlickerShiftyBar( RunIDs(r), w, ifi, Tex, NumMovieFrames, FlickerFrameIndices, BackgroundColor, ...
                        HorzBarHeightPix, VertBarWidthPix, HorzBarYsteps, VertBarXsteps, ...
                        HorFrameIx, HorFrameIxRev, VerFrameIx, VerFrameIxRev, RandShiftFrameIndices, Options );
            elseif strcmpi( Type, 'MovGrat' ) == 1
                [LoopTimestampHor,LoopTimestampVer] = ...
                    MovGratBar( RunIDs(r), w, ifi, Tex, NumMovieFrames, AngleFrameIndices, BackgroundColor, ...
                        HorzBarHeightPix, VertBarWidthPix, HorzBarYsteps, VertBarXsteps, ...
                        HorFrameIx, HorFrameIxRev, VerFrameIx, VerFrameIxRev, ...
                        xSteps, ySteps, MovieFrameIndices, Options );
            end
            StimSignal( Settings, Settings.test, [], 0 )
            Screen('FillRect',w, BackgroundColor);
            Screen('Flip', w);
            Screen('FillRect',w, BackgroundColor);
            TimeFinished = toc;
            disp(['  * Stimulus finished (t=' num2str(TimeFinished) ')']);

        end
%         % Back to full frame
%         set(Settings.CamId, 'ROIPosition', Settings.ROIorig );

        %% properly close the textures and screen
        disp('Closing screenbuffer...');    
        LoadIdentityClut(w);
        Screen('CloseAll');
        Priority(0);
        ShowCursor;
            
    catch
        LoadIdentityClut(w);
        Screen('CloseAll');
        Priority(0);
        ShowCursor;
        clearvars -global -except Ard
        clearvars;
        % We throw the error again so the user sees the error description.
        psychrethrow(psychlasterror);
        % Bye bye!
        return;
    end
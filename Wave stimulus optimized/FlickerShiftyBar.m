function [LoopTimestampHor,LoopTimestampVer] = ...
    FlickerShiftyBar( Direction, w, ifi, Tex, NumMovieFrames, FlickerFrameIndices, BackgroundColor, ...
                HorzBarHeightPix, VertBarWidthPix, HorzBarYsteps, VertBarXsteps, ...
                HorFrameIx, HorFrameIxRev, VerFrameIx, VerFrameIxRev, RandShiftFrameIndices, Options )
    
    % Flip a grey screen to get time of retrace
    Screen('FillRect', w, BackgroundColor);
    vbl = Screen('Flip', w );

    % set source coordinates Horizontal bar
    if Direction == 1 || Direction == 2 || Direction == 5 || Direction == 6
        hx1 = 1; hy1 = 0;
        hx2 = Options.scrx;
        hy2 = HorzBarHeightPix;
    end
    
    % set source coordinates Vertical bar
    if Direction == 3 || Direction == 4 || Direction == 5 || Direction == 6
        vx1 = 1; vy1 = 0;
        vx2 = VertBarWidthPix;
        vy2 = Options.scry;
    end
    
    % enable random shifts in black/white flicker patches
    RandVerShift = round( HorzBarHeightPix .* rand(1,NumMovieFrames) );
    RandHorShift = round( VertBarWidthPix .* rand(1,NumMovieFrames) );
    
    % display movie in coordinate window
    i = 0;
    LoopCountHor = 0;
    LoopCountVer = 0;
    LoopTimestampHor = [];
    LoopTimestampVer = [];
    if Direction == 1 % Horizontal down
        while i < NumMovieFrames
            i = i + 1;
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [hx1+RandHorShift(RandShiftFrameIndices(i)) hy1 + HorzBarYsteps(HorFrameIx(i))+RandVerShift(RandShiftFrameIndices(i)) ...
                 hx2+RandHorShift(RandShiftFrameIndices(i)) hy2 + HorzBarYsteps(HorFrameIx(i))+RandVerShift(RandShiftFrameIndices(i))],...
                 [hx1 HorzBarYsteps(HorFrameIx(i)) hx2 HorzBarYsteps(HorFrameIx(i)) + HorzBarHeightPix ] );
            vbl = Screen('Flip', w, vbl + 0.5 * ifi);
            if HorFrameIx(i) == max(HorFrameIx)
                Screen('FillRect', w, BackgroundColor);
                vbl = Screen('Flip', w );
                pause(Options.RepWait);
                LoopCountHor = LoopCountHor + 1;
                LoopTimestampHor(LoopCountHor) = toc;
                disp(['Loop ' num2str(LoopCountHor) ' (horizontal), time = ' num2str(LoopTimestampHor(LoopCountHor))]);
            end
            Exit_if_Shift_LOCAL
        end
    elseif Direction == 2 % Horizontal up
        while i < NumMovieFrames
            i = i + 1;
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [hx1+RandHorShift(RandShiftFrameIndices(i)) hy1 + HorzBarYsteps(HorFrameIxRev(i))+RandVerShift(RandShiftFrameIndices(i)) ...
                 hx2+RandHorShift(RandShiftFrameIndices(i)) hy2 + HorzBarYsteps(HorFrameIxRev(i))+RandVerShift(RandShiftFrameIndices(i))],...
                 [hx1 HorzBarYsteps(HorFrameIxRev(i)) hx2 HorzBarYsteps(HorFrameIxRev(i)) + HorzBarHeightPix ] );
            vbl = Screen('Flip', w, vbl + 0.5 * ifi);
            if HorFrameIx(i) == max(HorFrameIx)
                Screen('FillRect', w, BackgroundColor);
                vbl = Screen('Flip', w );
                pause(Options.RepWait);
                LoopCountHor = LoopCountHor + 1;
                LoopTimestampHor(LoopCountHor) = toc;
                disp(['Loop ' num2str(LoopCountHor) ' (horizontal), time = ' num2str(LoopTimestampHor(LoopCountHor))]);
            end
            Exit_if_Shift_LOCAL
        end
    elseif Direction == 3 % Vertical to right of screen
        while i < NumMovieFrames
            i = i + 1;
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [vx1 + VertBarXsteps(VerFrameIx(i))+RandHorShift(RandShiftFrameIndices(i)) vy1+RandVerShift(RandShiftFrameIndices(i)) ...
                 vx2 + VertBarXsteps(VerFrameIx(i))+RandHorShift(RandShiftFrameIndices(i)) vy2+RandVerShift(RandShiftFrameIndices(i)) ],...
                 [VertBarXsteps(VerFrameIx(i)) vy1 VertBarXsteps(VerFrameIx(i)) + VertBarWidthPix vy2 ] );
            vbl = Screen('Flip', w, vbl + 0.5 * ifi);
            if VerFrameIx(i) == max(VerFrameIx)
                Screen('FillRect', w, BackgroundColor);
                vbl = Screen('Flip', w );
                pause(Options.RepWait);
                LoopCountVer = LoopCountVer + 1;
                LoopTimestampVer(LoopCountVer) = toc;
                disp(['Loop ' num2str(LoopCountVer) ' (vertical), time = ' num2str(LoopTimestampVer(LoopCountVer))]);
            end
            Exit_if_Shift_LOCAL
        end
    elseif Direction == 4 % Vertical to left of screen
        while i < NumMovieFrames
            i = i + 1;
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [vx1 + VertBarXsteps(VerFrameIxRev(i))+RandHorShift(RandShiftFrameIndices(i)) vy1+RandVerShift(RandShiftFrameIndices(i)) ...
                 vx2 + VertBarXsteps(VerFrameIxRev(i))+RandHorShift(RandShiftFrameIndices(i)) vy2+RandVerShift(RandShiftFrameIndices(i)) ],...
                 [VertBarXsteps(VerFrameIxRev(i)) vy1 VertBarXsteps(VerFrameIxRev(i)) + VertBarWidthPix vy2 ] );
            vbl = Screen('Flip', w, vbl + 0.5 * ifi);
            if VerFrameIx(i) == max(VerFrameIx)
                Screen('FillRect', w, BackgroundColor);
                vbl = Screen('Flip', w );
                pause(Options.RepWait);
                LoopCountVer = LoopCountVer + 1;
                LoopTimestampVer(LoopCountVer) = toc;
                disp(['Loop ' num2str(LoopCountVer) ' (vertical), time = ' num2str(LoopTimestampVer(LoopCountVer))]);
            end
            Exit_if_Shift_LOCAL
        end
    elseif Direction == 5 % Horizontal down & Vertical to right of screen
        while i < NumMovieFrames
            i = i + 1;
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [hx1+RandHorShift(RandShiftFrameIndices(i)) hy1 + HorzBarYsteps(HorFrameIx(i))+RandVerShift(RandShiftFrameIndices(i)) ...
                 hx2+RandHorShift(RandShiftFrameIndices(i)) hy2 + HorzBarYsteps(HorFrameIx(i))+RandVerShift(RandShiftFrameIndices(i))],...
                 [hx1 HorzBarYsteps(HorFrameIx(i)) hx2 HorzBarYsteps(HorFrameIx(i)) + HorzBarHeightPix ] );
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [vx1 + VertBarXsteps(VerFrameIx(i))+RandHorShift(RandShiftFrameIndices(i)) vy1+RandVerShift(RandShiftFrameIndices(i)) ...
                 vx2 + VertBarXsteps(VerFrameIx(i))+RandHorShift(RandShiftFrameIndices(i)) vy2+RandVerShift(RandShiftFrameIndices(i)) ],...
                 [VertBarXsteps(VerFrameIx(i)) vy1 VertBarXsteps(VerFrameIx(i)) + VertBarWidthPix vy2 ] );
            vbl = Screen('Flip', w, vbl + 0.5 * ifi);
            if HorFrameIx(i) == 1
                LoopCountHor = LoopCountHor + 1;
                LoopTimestampHor(LoopCountHor) = toc;
                disp(['Loop ' num2str(LoopCountHor) ' (horizontal), time = ' num2str(LoopTimestampHor(LoopCountHor))]);
            end
            if VerFrameIx(i) == 1
                LoopCountVer = LoopCountVer + 1;
                LoopTimestampVer(LoopCountVer) = toc;
                disp(['Loop ' num2str(LoopCountVer) ' (vertical), time = ' num2str(LoopTimestampVer(LoopCountVer))]);
            end
            Exit_if_Shift_LOCAL
        end
    elseif Direction == 6
        while i < NumMovieFrames
            i = i + 1;
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [hx1+RandHorShift(RandShiftFrameIndices(i)) hy1 + HorzBarYsteps(HorFrameIxRev(i))+RandVerShift(RandShiftFrameIndices(i)) ...
                 hx2+RandHorShift(RandShiftFrameIndices(i)) hy2 + HorzBarYsteps(HorFrameIxRev(i))+RandVerShift(RandShiftFrameIndices(i))],...
                 [hx1 HorzBarYsteps(HorFrameIxRev(i)) hx2 HorzBarYsteps(HorFrameIxRev(i)) + HorzBarHeightPix ] );
            Screen('DrawTexture', w, Tex(FlickerFrameIndices(i)), ...
                [vx1 + VertBarXsteps(VerFrameIxRev(i))+RandHorShift(RandShiftFrameIndices(i)) vy1+RandVerShift(RandShiftFrameIndices(i)) ...
                vx2 + VertBarXsteps(VerFrameIxRev(i))+RandHorShift(RandShiftFrameIndices(i)) vy2+RandVerShift(RandShiftFrameIndices(i)) ],...
                 [VertBarXsteps(VerFrameIxRev(i)) vy1 VertBarXsteps(VerFrameIxRev(i)) + VertBarWidthPix vy2 ] );
            vbl = Screen('Flip', w, vbl + 0.5 * ifi);
            if HorFrameIx(i) == 1
                LoopCountHor = LoopCountHor + 1;
                LoopTimestampHor(LoopCountHor) = toc;
                disp(['Loop ' num2str(LoopCountHor) ' (horizontal), time = ' num2str(LoopTimestampHor(LoopCountHor))]);
            end
            if VerFrameIx(i) == 1
                LoopCountVer = LoopCountVer + 1;
                LoopTimestampVer(LoopCountVer) = toc;
                disp(['Loop ' num2str(LoopCountVer) ' (vertical), time = ' num2str(LoopTimestampVer(LoopCountVer))]);
            end
            Exit_if_Shift_LOCAL
        end
    end
    
    % Flip a grey screen and be done with it!
    Screen('FillRect', w, BackgroundColor);
    vbl = Screen('Flip', w );
    
end



function Exit_if_Shift_LOCAL

global Ard

if     strfind(eval('computer'),'WIN')
    yesKey = KbName('shift');
elseif strfind(eval('computer'),'LNX')
    yesKey = KbName('Shift_L');
elseif strfind(eval('computer'),'MAC')
    yesKey = KbName('LeftShift');
end
[~,~,keyCode] = KbCheck;

    if keyCode(yesKey)
        if ~isempty(Ard)
            try
                ArdPinNr = evalin('caller', 'Param.ArdPinNr');
            catch
                ArdPinNr = evalin('caller', 'Options.ArdPinNr');
            end
            Ard.digitalWrite(ArdPinNr,0); % set 0V
        end
        Screen('CloseAll');
        Priority(0);
        ShowCursor;
        return
    end
    
end
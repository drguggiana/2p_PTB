function [RBS, vbl] = PresentRandomBarStereogram( nidaq, screen, vbl, Param, RBS ) %#ok<INUSL>

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));

win=screen(1);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
Screen('SelectStereoDrawBuffer', win, 0);
Screen('FillRect', win, RBS.BackgroundLuminance);
Screen('SelectStereoDrawBuffer', win, 1);
Screen('FillRect', win, RBS.BackgroundLuminance);
vbl = Screen('Flip', win);

test=Param.test;
tempo=[];
cnt = RBS.cnt;
stimseq      = RBS.stimseq';
seqangles    = RBS.seqangles';
seqIntervals = RBS.seqIntervals';
angle_cart = seqangles(stimseq(cnt));
angle0 = angle_cart+90;
if angle0 >= 180
    angle1 = angle0 - 180 + RBS.offset_rot_deg(1);
    angle2 = angle0 - 180 + RBS.offset_rot_deg(2);
else
    angle1 = angle0 + RBS.offset_rot_deg(1);
    angle2 = angle0 + RBS.offset_rot_deg(2);
end

% angle = seqangles(RBS.cnt)
% -angle+270
    m = tan((-angle0+270)*pi/180);
    C = nchoosek([1:RBS.nBars],2);

disp([' --> Angle (cart) = ' num2str(mod(angle_cart,180)) ]);

list_intervals = 1 : RBS.n_disparity_intervals;
whichInterval = list_intervals(seqIntervals(stimseq(cnt)));

bin_ranges_deg = RBS.bin_ranges_deg(whichInterval,:);

%% generate and start stimulus
ticstart = tic;
vblendtime = vbl + RBS.stimulus_time - Param.ifi;
frameNr = 0;
patternNr = 0;

while vbl < vblendtime
    
    frameNr = frameNr + 1;
    
    vbl_patternend = vbl + RBS.pattern_time - Param.ifi;
    if vbl_patternend > vblendtime
        vbl_patternend = vblendtime;
    end
    go_pattern = 1;
    patternNr = patternNr + 1;
    
%     dimX = round( BM.barWidth_px  / 2 );
    dimY = round( RBS.barHeight_px / 2 );
%     baseRect = [-dimX -dimY, dimX dimY];
    
    % get nBars random numbers between margin and
    % Param.screenRes(1)-margin:
    repeat = 1;
    margin = max(RBS.barWidth_px);
%     while repeat==1
%         posXs = ((Param.screenRes(1)-margin)-margin)*rand(BM.nBars,1) + margin;
%         posYs = ((Param.screenRes(2)-margin)-margin)*rand(BM.nBars,1) + margin;
%         if all(abs(diff(posXs)) > BM.barWidth_px*3 ) &&...
%             all(abs(diff(posYs)) > BM.barWidth_px*3 )
%         posXs
%         posYs
%         abs(diff(posXs))
%         abs(diff(posYs))
%             repeat=0; break %#ok<NASGU>
%         end
%     end

    while repeat==1
        d=zeros(RBS.nBars,1);
        posXs = round(((Param.screenRes(1)-margin)-margin)*rand(RBS.nBars,1) + margin);
        posYs = round(((Param.screenRes(2)-margin)-margin)*rand(RBS.nBars,1) + margin);
%         posYs = Param.screenRes(2)-posYs;
        if     angle0 == 90 || angle0 == 270
            if all(abs(diff(posYs)) > max(RBS.barWidth_px)*2 )
                repeat=0; break
            end
        elseif angle0 == 0 || angle0 == 180
            if all(abs(diff(posXs)) > max(RBS.barWidth_px)*2 )
                repeat=0; break
            end
        else
            for c = 1:size(C,1)
                i1 = C(c,1); i2 = C(c,2);
    %             d(c) = abs(posYs(i2)-posYs(i1)-m*(posXs(i2)-posXs(i1))) / ...
    %                     sqrt(m^2+1) ;
                d(c) = abs(-posYs(i2)+posYs(i1)-m*(posXs(i2)-posXs(i1))) / ...
                        sqrt(m^2+1) ;
            end
            if all(abs(diff(d)) > max(RBS.barWidth_px)*2 )
                repeat=0; break %#ok<NASGU>
            end
        end
    end
%     posYs = ones(1,BM.nBars) .* (Param.screenRes(2) / 2);
%     posYs = Param.screenRes(2)-posYs;

    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 0);

    for n = 1 : RBS.nBars
        dimX(n) = round((RBS.barWidth_px(2)-RBS.barWidth_px(1))*rand(1,1) + RBS.barWidth_px(1));
        baseRect(n,:) = [-dimX(n) -dimY, dimX(n) dimY];
        % Get the current squares position ans rotation angle
        posX = posXs(n);
        posY = posYs(n);
        % Translate, rotate, re-tranlate and then draw our square
        Screen('glPushMatrix', win)
        Screen('glTranslate', win, posX, posY)
        Screen('glRotate', win, angle1, 0, 0);
        Screen('glTranslate', win, -posX, -posY)
        Screen('FillRect', win, RBS.color_list(:,n),...
            CenterRectOnPoint(baseRect(n,:), posX, posY));
        Screen('glPopMatrix', win)
        
    end

    
    % Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', win, 1);
    
    
    disparity_deg = (bin_ranges_deg(2)-bin_ranges_deg(1))*rand(1) + bin_ranges_deg(1);
    dx =  disparity_deg*sin(angle_cart*pi/180);
    dy = -disparity_deg*cos(angle_cart*pi/180);

    
    disp([ ' ----> Range disparity [' num2str(bin_ranges_deg) '] - value: ' num2str(disparity_deg,'%2.1f')]);
    
   


        for n = 1 : RBS.nBars
            % Get the current squares position ans rotation angle
            posX = posXs(n) + dx;
            posY = posYs(n) + dy;
            % Translate, rotate, re-tranlate and then draw our square
            Screen('glPushMatrix', win)
            Screen('glTranslate', win, posX, posY)
            Screen('glRotate', win, angle2, 0, 0);
            Screen('glTranslate', win, -posX, -posY)
            Screen('FillRect', win, RBS.color_list(:,n),...
                CenterRectOnPoint(baseRect(n,:), posX, posY));
            Screen('glPopMatrix', win)
        end

        dontclear = 1;
        while go_pattern == 1
            if frameNr == 1
                StimSignal( Param, test, nidaq, RBS.stim_id );
            end
            vbl = Screen('Flip', win, vbl + 0.5*Param.ifi, dontclear);
            if vbl > vbl_patternend
                dontclear = 0;
                vbl = Screen('Flip', win, vbl + 0.5*Param.ifi, dontclear);
                go_pattern = 0;
                tempo(patternNr) = toc(ticstart); %#ok<AGROW>
            end
            Exit_If_Shift(test,nidaq);
        end
        

    
end

StimSignal( Param, test, nidaq, 0 );
% end stimulus


%% post-stimulus
ticpoststim = tic;
Screen('FillRect', win, RBS.BackgroundLuminance);

for i = 1 : RBS.frames_poststim
    vbl = Screen('Flip', win);
    Exit_If_Shift(test,nidaq);
end
tempo(end+1) = toc(ticpoststim);
% end post-stimulus
RBS.tempo{cnt}=tempo;


if exist('scriptName','var') && RBS.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
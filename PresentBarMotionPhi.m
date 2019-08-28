function [BM, vbl] = PresentBarMotionPhi( nidaq, screen, vbl, Param, BM, control )

[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
if ~exist('control','var') || isempty(control)
    control = 0;
end
win=screen(1);
Screen('FillRect', win, BM.BackgroundLuminance);
vbl = Screen('Flip', win);
Screen('BlendFunction', win, GL_ONE, GL_ZERO);
test=Param.test;
tempo=[0,0];
seqangles = BM.seqangles';
angle = seqangles(BM.cnt)
% -angle+270
    m = tan((-angle+270)*pi/180);
    C = nchoosek([1:BM.nBars],2);

ticstart = tic;
vblendtime = vbl + BM.stimulus_time - Param.ifi;
go_on=1;

while go_on == 1  
    
    
    
%     dimX = round( BM.barWidth_px  / 2 );
    dimY = round( BM.barHeight_px / 2 );
%     baseRect = [-dimX -dimY, dimX dimY];
    
    % get nBars random numbers between margin and
    % Param.screenRes(1)-margin:
    repeat = 1;
    margin = max(BM.barWidth_px);
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
        d=zeros(BM.nBars,1);
        posXs = round(((Param.screenRes(1)-margin)-margin)*rand(BM.nBars,1) + margin);
        posYs = round(((Param.screenRes(2)-margin)-margin)*rand(BM.nBars,1) + margin);
%         posYs = Param.screenRes(2)-posYs;
        if     angle == 90 || angle == 270
            if all(abs(diff(posYs)) > max(BM.barWidth_px)*2 )
                repeat=0; break
            end
        elseif angle == 0 || angle == 180
            if all(abs(diff(posXs)) > max(BM.barWidth_px)*2 )
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
            if all(abs(diff(d)) > max(BM.barWidth_px)*2 )
                repeat=0; break %#ok<NASGU>
            end
        end
    end
%     posYs = ones(1,BM.nBars) .* (Param.screenRes(2) / 2);
%     posYs = Param.screenRes(2)-posYs;
    
    for n = 1 : BM.nBars
        dimX(n) = round((BM.barWidth_px(2)-BM.barWidth_px(1))*rand(1,1) + BM.barWidth_px(1));
        baseRect(n,:) = [-dimX(n) -dimY, dimX(n) dimY];
        % Get the current squares position ans rotation angle
        posX = posXs(n);
        posY = posYs(n);
        % Translate, rotate, re-tranlate and then draw our square
        Screen('glPushMatrix', win)
        Screen('glTranslate', win, posX, posY)
        Screen('glRotate', win, angle, 0, 0);
        Screen('glTranslate', win, -posX, -posY)
        Screen('FillRect', win, BM.color_list(:,n),...
            CenterRectOnPoint(baseRect(n,:), posX, posY));
        Screen('glPopMatrix', win)
        
%         Screen('DrawDots',win, [posX;posY], 50, 0, [0 0], 1);
    end
    
    tic
    StimSignal( Param, test, nidaq, stimtype.stim_id )
    dontclear=1;
    for ii = 1:BM.phi_interval_fr-1
		if Param.usePhotodiode 
			Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
		end
        vbl = Screen('Flip', win, [], dontclear);
%     figure;plot(posXs,-posYs,'o')
        if vbl > vblendtime
            go_on=0;
            break
        end
        Exit_If_Shift(test,nidaq);
    end
    vbl = Screen('Flip', win, [], 0);
    toc
    
    if control==0
        
        phi_jump = (BM.phi_jump(2)-BM.phi_jump(1))*rand(1) + BM.phi_jump(1);
    %     dx =  phi_jump * Param.pixperdeg;
    %     dy = 0;
    %     dy = -phi_jump * Param.pixperdeg;
        dx =  phi_jump*cos(angle*pi/180) * Param.pixperdeg;
    %     dy = -phi_jump*cos(angle*pi/180) * Param.pixperdeg;
        dy = +phi_jump*sin(angle*pi/180) * Param.pixperdeg;

        for n = 1 : BM.nBars
            % Get the current squares position ans rotation angle
            posX = posXs(n) + dx;
            posY = posYs(n) + dy;
            % Translate, rotate, re-tranlate and then draw our square
            Screen('glPushMatrix', win)
            Screen('glTranslate', win, posX, posY)
            Screen('glRotate', win, angle, 0, 0);
            Screen('glTranslate', win, -posX, -posY)
            Screen('FillRect', win, BM.color_list(:,n),...
                CenterRectOnPoint(baseRect(n,:), posX, posY));
            Screen('glPopMatrix', win)
        end

        tic
        StimSignal( Param, test, nidaq, BM.stim_id )
        dontclear=1;
        for ii = 1:BM.phi_interval_fr-1
			if Param.usePhotodiode 
				Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
			end
            vbl = Screen('Flip', win, [], dontclear);
            if vbl > vblendtime
                go_on=0;
                break
            end
            Exit_If_Shift(test,nidaq);
        end
        vbl = Screen('Flip', win, [], 0);
        toc
        
    elseif control==1
        
        tic
        Screen('FillRect', win, BM.BackgroundLuminance);
        for ii = 1:BM.phi_interval_fr

            vbl = Screen('Flip', win);
            if vbl > vblendtime
                go_on=0;
                break
            end
            Exit_If_Shift(test,nidaq);
        end
        toc
        
    end
    
end

tempo(1) = toc(ticstart);
% end stimulus

ticpoststim = tic;
StimSignal( Param, test, nidaq, 0 )

% post-stimulus
Screen('FillRect', win, BM.BackgroundLuminance);

for n = 1:BM.frames_poststim
    vbl = Screen('Flip', win);
    Exit_If_Shift(test,nidaq);
end
tempo(2) = toc(ticpoststim);
% end post-stimulus
BM.tempo{BM.cnt}=tempo;

if exist('scriptName','var') && BM.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end
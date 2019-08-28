
function [stimtype, vbl] = PresentPhiMotionGrating( nidaq, screen, vbl, Param, stimtype )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
tic;

if Param.stereoMode == 4 || Param.stereoMode == 1
    win=screen(1);
    Screen('BlendFunction', win, GL_ONE, GL_ZERO);
    Screen('SelectStereoDrawBuffer', win, 0);
    Screen('FillRect', win, stimtype.BackgroundLuminance);
    Screen('SelectStereoDrawBuffer', win, 1);
    Screen('FillRect', win, stimtype.BackgroundLuminance);
    vbl = Screen('Flip', win);
elseif Param.stereoMode == 0
    for i=1:length(screen)
        Screen('BlendFunction', screen(i), GL_ONE, GL_ZERO);
    end
    if     Param.Dichoptic
        Screen('FillRect', screen(2), stimtype.BackgroundLuminance);
        vbl = Screen('Flip', screen(2));
    end
end

test=Param.test;
tempo=[];
cnt=stimtype.cnt;
stimseq=stimtype.stimseq; % (check that stimseq columns => reps)
seqangles = stimtype.seqangles; % (check that seqangles columns => different orientations)

angle0 = seqangles(stimseq(cnt));
anglecart = mod(angle0-90, 180);
if angle0 >= 180
    angle1 = angle0 - 180 + stimtype.offset_rot_deg(1);
    angle2 = angle0 - 180 + stimtype.offset_rot_deg(2);
%     phaseincrement = - stimtype.phaseincrement;
else
    angle1 = angle0 + stimtype.offset_rot_deg(1);
    angle2 = angle0 + stimtype.offset_rot_deg(2);
%     phaseincrement = + stimtype.phaseincrement;
end

phase1 = stimtype.seqphases(stimseq(cnt),1);
phase2 = stimtype.seqphases(stimseq(cnt),2);
fprintf('  * Orientation=%3.0f, (cartesian %3.0f); Phase1=%3.0f, Phase2=%3.0f \n',angle0,anglecart,phase1,phase2);
phase1 = phase1 + abs(cos(rad(angle1)))*stimtype.offset_phase1_deg(1) + abs(sin(rad(angle1)))*stimtype.offset_phase1_deg(2);
phase2 = phase2 + abs(cos(rad(angle2)))*stimtype.offset_phase2_deg(1) + abs(sin(rad(angle2)))*stimtype.offset_phase2_deg(2);


vblendtime = vbl + stimtype.stimulus_time - Param.ifi;  %Ale: I put -ifi to obtain an exact stimulus time.
go_on = 1;
frameNr=0;

while go_on == 1  

% %     for p = repmat([1 2],1,20)
        if go_on==1
            tic
% %             if     p == 1
% %                 % to bring the white bar to the center of the RF that you
% %                 % set with PMG.offset_phase1_perc:
% %                     phase1 = stimtype.phase1(seqphases(stimseq(cnt))) + abs(cos(rad(angle1)))*stimtype.offset_phase1_deg(1) + abs(sin(rad(angle1)))*stimtype.offset_phase1_deg(2)
% %                 if Param.Dichoptic
% %                     phase2 = stimtype.phase2(seqphases(stimseq(cnt))) + abs(cos(rad(angle2)))*stimtype.offset_phase2_deg(1) + abs(sin(rad(angle2)))*stimtype.offset_phase2_deg(2)
% %                 end
% %             elseif p == 2
% %                     phase1 = stimtype.phase1(seqphases(stimseq(cnt))) + abs(cos(rad(angle1)))*stimtype.offset_phase1_deg(1) + abs(sin(rad(angle1)))*stimtype.offset_phase1_deg(2) +180;
% %                 if Param.Dichoptic
% %                     phase2 = stimtype.phase2(seqphases(stimseq(cnt))) + abs(cos(rad(angle2)))*stimtype.offset_phase2_deg(1) + abs(sin(rad(angle2)))*stimtype.offset_phase2_deg(2) +180;
% %                 end
% %             end

% %             vblendpatch = vbl + stimtype.patch_time - Param.ifi;
            while vbl < vblendtime
                
                
                amplitude1 = stimtype.amplitude * cos(2*pi*stimtype.cyclespersecond*frameNr*Param.ifi);
                % amplitude2 is shifted half period compared to amplitude1:
                amplitude2 = stimtype.amplitude * sin(2*pi*stimtype.cyclespersecond*frameNr*Param.ifi);
                frameNr = frameNr+1;
                if frameNr == 1
                    StimSignal( Param, test, nidaq, stimtype.stim_id )
                end
                if     Param.stereoMode == 4 || Param.stereoMode == 1
                        % Select left-eye image buffer for drawing:
                        Screen('SelectStereoDrawBuffer', win, 0);
                        Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingsize, angle1, [], [], [], [], Param.rotateMode, [phase1, stimtype.freq_cppx(1), amplitude1, 0]);
                        if Param.usePhotodiode
                            Screen('SelectStereoDrawBuffer', win, Param.Photodiode.screen); 
                            Screen('FillRect', win , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
                        end
                        Screen('SelectStereoDrawBuffer', win, 1);
                        Screen('DrawTexture', win, stimtype.gratingtex, [], stimtype.gratingsize, angle2, [], [], [], [], Param.rotateMode, [phase2, stimtype.freq_cppx(2), amplitude2, 0]);
                        Screen('DrawingFinished', win);
                        vbl = Screen('Flip', win, vbl + 0.5 * Param.ifi);
                elseif Param.stereoMode == 0
                    if Param.Dichoptic
                        Screen('DrawTexture', screen(1), stimtype.gratingtex, [], stimtype.gratingsize, angle1, [], [], [], [], Param.rotateMode, [phase1, stimtype.freq_cppx(1), amplitude1, 0]);
                        Screen('DrawTexture', screen(2), stimtype.gratingtex, [], stimtype.gratingsize, angle2, [], [], [], [], Param.rotateMode, [phase2, stimtype.freq_cppx(2), amplitude2, 0]);
                        if Param.usePhotodiode
                            Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
                        end
    % %                     Screen('FillRect', screen(2), stimtype.BackgroundLuminance );
    %                     vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
    %                     vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi );
                        vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
                        vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi );
    %                     vbl=vbl2;
                    else
    %                     Screen('DrawTexture', screen(1), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase1, stimtype.freq_cppx(1), amplitude, 0]);
    % 					if Param.usePhotodiode
    % 						Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
    % 					end
    %                     vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
                    end
                end
                

                if vbl > vblendtime
                    go_on=0;
                    break
                end
                Exit_If_Shift(test,nidaq);
            end
            tempo = [tempo toc]; %#ok<*AGROW>

             % interpatch #1
% %                 if stimtype.interpatch_time > 0 && go_on==1
% %                     tic
% %                 % StimSignal( Param, test, nidaq, 0 )
% %                 % # frames of blank screen between first and second flashed
% %                 % gratings (don't know why I need -1 to get the precise
% %                 % timing):
% % %                     for ff=1:stimtype.frame_interpatch-1    
% % %                     while toc < stimtype.interpatch_time
% %                     vblendinterpatch = vbl + stimtype.interpatch_time - Param.ifi;
% %                     while vbl < vblendinterpatch
% % %                         for i=1:length(screen)
% %                         if Param.Dichoptic
% %                             Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
% % 							if Param.usePhotodiode
% % 								Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
% % 							end
% %                             Screen('FillRect', screen(2), stimtype.BackgroundLuminance );
% % %                             vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi);
% % %                             vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi);
% %                             vbl  = Screen('Flip', screen(1) );
% %                             vbl2 = Screen('Flip', screen(2) );
% %                             vbl=vbl2;
% %                         else
% %                             Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
% % 							if Param.usePhotodiode
% % 								Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
% % 							end
% %                             vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi);
% %                         end
% %                             if vbl > vblendtime
% %                                 go_on=0;
% %                                 break
% %                             end
% %                     end
% %                     tempo = [tempo toc];
% %                 end

            % post-patch blank screen for screen1 and stim for screen2    
% %             if go_on==1
% %                 tic
% % %                 for f = 1 : stimtype.patch_time_fr
% % %                 while toc < stimtype.patch_time
% %                 vblendpatch = vbl + stimtype.patch_time - Param.ifi;
% %                 while vbl < vblendpatch
% %                     if Param.Dichoptic
% %                         Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
% % 						if Param.usePhotodiode
% % 							Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
% % 						end
% %                         Screen('DrawTexture', screen(2), stimtype.gratingtex, [], stimtype.gratingsize, angle+stimtype.offset_rot_deg, [], [], [], [], Param.rotateMode, [phase2, stimtype.freq_cppx(2), stimtype.amplitude, 0]);
% % %                         vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi);
% % %                         vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi);
% %                         vbl  = Screen('Flip', screen(1) );
% %                         vbl2 = Screen('Flip', screen(2) );
% %                         vbl=vbl2;
% %                     else
% %                         Screen('DrawTexture', screen(1), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase1+phaseincrement, stimtype.freq_cppx(1), stimtype.amplitude, 0]);
% % 						if Param.usePhotodiode
% % 							Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
% % 						end
% %                         vbl = Screen('Flip', screen(1));
% %                     end
% % 
% %                     if vbl > vblendtime
% %                         go_on=0;
% %                         break
% %                     end
% %                     Exit_If_Shift(test,nidaq);
% %                 end
% %                 tempo = [tempo toc];
% %             end
            
            % interpatch #2
% %                 if stimtype.interpatch_time > 0 && go_on==1
% %                     tic
% %                 % StimSignal( Param, test, nidaq, 0 )
% %                 % # frames of blank screen between second and first flashed
% %                 % gratings (don't know why I need -1 to get the precise
% %                 % timing):
% % %                     for ff=1:stimtype.frame_interpatch-1   
% % %                     while toc < stimtype.interpatch_time
% %                     vblendinterpatch = vbl + stimtype.interpatch_time - Param.ifi;
% %                     while vbl < vblendinterpatch
% % %                         for i=1:length(screen)   
% %                             if Param.Dichoptic
% %                                 Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
% %                                 Screen('FillRect', screen(2), stimtype.BackgroundLuminance );
% % %                                 vbl = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
% % %                                 vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi );
% %                                 vbl  = Screen('Flip', screen(1) );
% %                                 vbl2 = Screen('Flip', screen(2) );
% %                                 vbl=vbl2;
% %                             else
% %                                 Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
% %                                 vbl = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
% %                             end
% %                             if vbl > vblendtime
% %                                 go_on=0;
% %                                 break
% %                             end
% % %                         end
% %                     end
% %                     tempo = [tempo toc];
% %                 end
                
        else
            break
        end
% %     end % for p = repmat([1 2],1,20)
end

tempo = [tempo toc];
StimSignal( Param, test, nidaq, 0 )

tic
for i=1:length(screen)
    Screen('FillRect', screen(i), stimtype.BackgroundLuminance);
end
if Param.stereoMode == 4 || Param.stereoMode == 1
    while  vbl < vblendtime + stimtype.poststim_time
        Exit_If_Shift(test,nidaq);
        vbl = Screen('Flip', win);  
    end
elseif Param.stereoMode == 0
    while  vbl < vblendtime + stimtype.poststim_time
        Exit_If_Shift(test,nidaq);
        vbl = Screen('Flip', screen(1));  
        if Param.Dichoptic
            Screen('Flip', screen(2));
        end
    end
end

tempo = [tempo toc];
stimtype.tempo{cnt}=tempo;


if exist('scriptName','var') && stimtype.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end

    function alpharad = rad(alphadeg)
        alpharad=alphadeg*pi/180;
    end
end
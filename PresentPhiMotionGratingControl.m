
function [stimtype, vbl] = PresentPhiMotionGratingControl( nidaq, screen, vbl, Param, stimtype, screen_id )
[scriptPath, scriptName] = fileparts(mfilename('fullpath'));
tic;
for i=1:length(screen)
    Screen('BlendFunction', screen(i), GL_ONE, GL_ZERO);
end
% if     Param.Dichoptic
%     Screen('FillRect', screen(2), stimtype.BackgroundLuminance);
%     vbl = Screen('Flip', screen(2));
% end
Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
Screen('Flip', screen(1));
Screen('FillRect', screen(2), stimtype.BackgroundLuminance );
vbl = Screen('Flip', screen(2));

test=Param.test;
tempo=[];
cnt=stimtype.cnt;
stimseq=stimtype.stimseq';
seqangles = stimtype.seqangles';
seqphases = stimtype.seqphases';
angle = seqangles(stimseq(cnt));
vblendtime = vbl + stimtype.stimulus_time - Param.ifi;  %Ale: I put -ifi to obtain an exact stimulus time, I dont know why exactly.
go_on = 1;

while go_on == 1  
    
    StimSignal( Param, test, nidaq, stimtype.stim_id )
    
    for p = repmat([1 2],1,40)
        if go_on==1
            tic
            if     p == 1
                if screen_id==1
                % to bring the white bar to the center of the RF that you
                % set with PMG.offset_phase1_perc:
                    phase = stimtype.phase1(seqphases(stimseq(cnt))) + abs(cos(rad(angle)))*stimtype.offset_phase1_deg(1) + abs(sin(rad(angle)))*stimtype.offset_phase1_deg(2);
                elseif screen_id==2
                    phase = stimtype.phase2(seqphases(stimseq(cnt))) + abs(cos(rad(angle)))*stimtype.offset_phase2_deg(1) + abs(sin(rad(angle)))*stimtype.offset_phase2_deg(2);
                end
            elseif p == 2
                if screen_id==1
                    phase = stimtype.phase1(seqphases(stimseq(cnt))) + abs(cos(rad(angle)))*stimtype.offset_phase1_deg(1) + abs(sin(rad(angle)))*stimtype.offset_phase1_deg(2) +180;
                elseif screen_id==2
                    phase = stimtype.phase2(seqphases(stimseq(cnt))) + abs(cos(rad(angle)))*stimtype.offset_phase2_deg(1) + abs(sin(rad(angle)))*stimtype.offset_phase2_deg(2) +180;
                end
            end
%             for f = 1 : stimtype.patch_time_fr-2
%             while toc < stimtype.patch_time
            vblendpatch = vbl + stimtype.patch_time - Param.ifi;
            while vbl < vblendpatch
                % Draw the grating, centered on the screen, with given rotation 'angle',
                % sine grating 'phase' shift and amplitude, rotating via set
                % 'rotateMode'. Note that we pad the last argument with a 4th
                % component, which is 0. This is required, as this argument must be a
                % vector with a number of components that is an integral multiple of 4,
                % i.e. in our case it must have 4 components:
%                 if Param.Dichoptic
%                     Screen('DrawTexture', screen(1), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase1, stimtype.freq_cppx(1), stimtype.amplitude, 0]);
%                     Screen('FillRect', screen(2), stimtype.BackgroundLuminance );
% %                     vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
% %                     vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi );
%                     vbl  = Screen('Flip', screen(1) );
%                     vbl2 = Screen('Flip', screen(2) );
%                     vbl=vbl2;
%                 else
                if screen_id==1
                    Screen('DrawTexture', screen(1), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase, stimtype.freq_cppx(1), stimtype.amplitude, 0]);
					if Param.usePhotodiode
						Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
					end
                    vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi );
                elseif screen_id==2
                    Screen('DrawTexture', screen(2), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase, stimtype.freq_cppx(2), stimtype.amplitude, 0]);
					if Param.usePhotodiode
						Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
					end
                    vbl  = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi );
                end
                if vbl > vblendtime
                    go_on=0;
                    break
                end
                Exit_If_Shift(test,nidaq);
            end
            tempo = [tempo toc]; %#ok<*AGROW>

             % interpatch #1
                if stimtype.interpatch_time > 0 && go_on==1
                    tic
%                 StimSignal( Param, test, nidaq, 0 )
                % # frames of blank screen between first and second flashed
                % gratings (don't know why I need -1 to get the precise
                % timing):
%                     for ff=1:stimtype.frame_interpatch-1    
%                     while toc < stimtype.interpatch_time
                    vblendinterpatch = vbl + stimtype.interpatch_time - Param.ifi;
                    while vbl < vblendinterpatch
%                         for i=1:length(screen)
                        if screen_id==1
                            Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
							if Param.usePhotodiode
								Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
							end
                            vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi);
                        elseif screen_id==2
                            Screen('FillRect', screen(2), stimtype.BackgroundLuminance );
							if Param.usePhotodiode
								Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
							end
                            vbl  = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi);
                        end
                            if vbl > vblendtime
                                go_on=0;
                                break
                            end
                    end
                    tempo = [tempo toc];
                end

            % post-patch blank screen for screen1 and stim for screen2    
            if go_on==1
                tic
%                 for f = 1 : stimtype.patch_time_fr
%                 while toc < stimtype.patch_time
                vblendpatch = vbl + stimtype.patch_time - Param.ifi;
                while vbl < vblendpatch
%                     if Param.Dichoptic
%                         Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
%                         Screen('DrawTexture', screen(2), stimtype.gratingtex, [], stimtype.gratingsize, angle+stimtype.offset_rot_deg, [], [], [], [], Param.rotateMode, [phase2, stimtype.freq_cppx(2), stimtype.amplitude, 0]);
% %                         vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi);
% %                         vbl2 = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi);
%                         vbl  = Screen('Flip', screen(1) );
%                         vbl2 = Screen('Flip', screen(2) );
%                         vbl=vbl2;
%                     else
                    if screen_id==1
                        Screen('DrawTexture', screen(1), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase+stimtype.phaseincrement, stimtype.freq_cppx(1), stimtype.amplitude, 0]);
						if Param.usePhotodiode
							Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
						end
                        vbl = Screen('Flip', screen(1));
                    elseif screen_id==2
                        Screen('DrawTexture', screen(2), stimtype.gratingtex, [], stimtype.gratingsize, angle, [], [], [], [], Param.rotateMode, [phase+stimtype.phaseincrement, stimtype.freq_cppx(2), stimtype.amplitude, 0]);
						if Param.usePhotodiode
							Screen('FillRect', screen(1) , Param.Photodiode.ColorOn, Param.Photodiode.Coord );
						end
                        vbl = Screen('Flip', screen(2));
                    end

                    if vbl > vblendtime
                        go_on=0;
                        break
                    end
                    Exit_If_Shift(test,nidaq);
                end
                tempo = [tempo toc];
            end
            
            % interpatch #2
                if stimtype.interpatch_time > 0 && go_on==1
                    tic
                    
%                 StimSignal( Param, test, nidaq, 0 )
                % # frames of blank screen between second and first flashed
                % gratings (don't know why I need -1 to get the precise
                % timing):
%                     for ff=1:stimtype.frame_interpatch-1   
%                     while toc < stimtype.interpatch_time
                    vblendinterpatch = vbl + stimtype.interpatch_time - Param.ifi;
                    while vbl < vblendinterpatch
%                         for i=1:length(screen)   
                        if screen_id==1
                            Screen('FillRect', screen(1), stimtype.BackgroundLuminance );
                            vbl  = Screen('Flip', screen(1), vbl + 0.5 * Param.ifi);
                        elseif screen_id==2
                            Screen('FillRect', screen(2), stimtype.BackgroundLuminance );
                            vbl  = Screen('Flip', screen(2), vbl + 0.5 * Param.ifi);
                        end
                            if vbl > vblendtime
                                go_on=0;
                                break
                            end
%                         end
                    end
                    tempo = [tempo toc];
                end
                
        else
            break
        end
    end
end
StimSignal( Param, test, nidaq, 0 )
tic
for i=1:length(screen)
Screen('FillRect', screen(i), stimtype.BackgroundLuminance);
end
while  vbl < vblendtime + stimtype.poststim_time
	Exit_If_Shift(test,nidaq);
    vbl = Screen('Flip', screen(1));  
    if Param.Dichoptic
        Screen('Flip', screen(2));
    end
end
tempo = [tempo toc];
stimtype.tempo{cnt}=tempo;


if exist('scriptName','var') && stimtype.cnt==1
    save_dir = evalin('caller', 'save_dir');
    BackupStimulusScript( scriptPath, scriptName, save_dir )
end

    function alpha2 = rad(alpha)
        alpha2=alpha*pi/180;
    end
end
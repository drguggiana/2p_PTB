function Exit_If_Shift(test,nidaq,only_shift)
% Exit presentation of stimuli when pressing Shift or closing laser
% shutter. StimFile is saved when exiting.
%
% When exiting, you get an error message like (just ignore it):
%
% Error in function DrawTexture: 	Invalid Window (or Texture) Index provided: It doesn't correspond to an open window or texture.
% Did you close it accidentally via Screen('Close') or Screen('CloseAll') ?
% Error using Screen
%
% The status of the shutter (open/close) is checked every
% 'IntervalShutterCheck' seconds, because querying inputSingleScan(nidaq)
% too often is demanding and causes screen tearing in the visual stimuli.

persistent TimerShutterCheck startTimer
global stim_fname
global Ard
global DG DG2 DG3 DG4 DG5 DG6
global DGD DGD2 DGD3
global DGF DGF2 DGF3
global DGID
global DGP
global SG
global RP RP2
global RPS RPS2
global RFM RFM2
global PMB
global PMG
global DGDphi
global PMGc1 PMGc2
global PG
global DM DM2 DM3
global RDS RDS2
global RBS RBS2
global NIS
global RDSB
global RDSE
global LO


IntervalShutterCheck = 0.2; % sec

if     strfind(eval('computer'),'WIN')
    yesKey = KbName('shift');
elseif strfind(eval('computer'),'LNX')
    yesKey = KbName('Shift_L');
elseif strfind(eval('computer'),'MAC')
    yesKey = KbName('LeftShift');
end
[~,~,keyCode] = KbCheck;

if nargin<3
    if keyCode(yesKey)
        saveStimFile
        if ~isempty(Ard)
            ArdPinNr = evalin('caller', 'Param.ArdPinNr');
            Ard.digitalWrite(ArdPinNr,0); % set 0V
        end
        ExitStim;
        return
    end
%     if  test == 0 && inputSingleScan(nidaq)<0.5;
%         saveStimFile
%         Screen('CloseAll');
%         Priority(0);
%         ShowCursor;
%         return
%     end
    if  test == 0
        if ~isempty(nidaq)
            if isempty(TimerShutterCheck)
                startTimer = tic;
            end
            TimerShutterCheck = toc(startTimer);
            % check shutter only after 'IntervalShutterCheck' sec, then reset
            % TimerShutterCheck:
            if TimerShutterCheck > IntervalShutterCheck
                if inputSingleScan(nidaq)<0.5;
                    saveStimFile
                    Screen('CloseAll');
                    Priority(0);
                    ShowCursor;
                    outputSingleScan(nidaq,0);
                    return
                end
                TimerShutterCheck = [];
            end
        end
        
    end
end
if nargin==3 && strcmp(only_shift,'only_shift')
    if keyCode(yesKey)
        ExitStim;
        return
    end
end


    function saveStimFile
        if ~isempty(DG)
            save(stim_fname, '-append','DG');
        end
        if ~isempty(LO)
            save(stim_fname, '-append','LO');
        end
        if ~isempty(DG2)
            save(stim_fname, '-append','DG2');
        end
        if ~isempty(DG3)
            save(stim_fname, '-append','DG3');
        end
        if ~isempty(DG4)
            save(stim_fname, '-append','DG4');
        end
        if ~isempty(DG5)
            save(stim_fname, '-append','DG5');
        end
        if ~isempty(DG6)
            save(stim_fname, '-append','DG6');
        end
        if ~isempty(DGD)
            save(stim_fname, '-append','DGD');
        end
        if ~isempty(DGD2)
            save(stim_fname, '-append','DGD2');
        end
        if ~isempty(DGD3)
            save(stim_fname, '-append','DGD3');
        end
        if ~isempty(DGF)
            save(stim_fname, '-append','DGF');
        end
        if ~isempty(DGF2)
            save(stim_fname, '-append','DGF2');
        end
        if ~isempty(DGF3)
            save(stim_fname, '-append','DGF3');
        end
        if ~isempty(DGID)
            save(stim_fname, '-append','DGID');
        end
        if ~isempty(DGP)
            save(stim_fname, '-append','DGP');
        end
        if ~isempty(SG)
            save(stim_fname, '-append','SG');
        end
        if ~isempty(RP)
            save(stim_fname, '-append','RP');
        end
        if ~isempty(RP2)
            save(stim_fname, '-append','RP2');
        end
        if ~isempty(RPS)
            save(stim_fname, '-append','RPS');
        end
        if ~isempty(RPS2)
            save(stim_fname, '-append','RPS2');
        end
        if ~isempty(RFM)
            save(stim_fname, '-append','RFM');
        end
        if ~isempty(RFM2)
            save(stim_fname, '-append','RFM2');
        end
        if ~isempty(PMG)
            save(stim_fname, '-append','PMG');
        end
        if ~isempty(DGDphi)
            save(stim_fname, '-append','DGDphi');
        end
        if ~isempty(PMGc1)
            save(stim_fname, '-append','PMGc1');
        end
        if ~isempty(PMGc2)
            save(stim_fname, '-append','PMGc2');
        end
        if ~isempty(PMB)
            save(stim_fname, '-append','PMB');
        end
        if ~isempty(PG)
            save(stim_fname, '-append','PG');
        end
        if ~isempty(DM)
            save(stim_fname, '-append','DM');
        end
        if ~isempty(DM2)
            save(stim_fname, '-append','DM2');
        end
        if ~isempty(DM3)
            save(stim_fname, '-append','DM3');
        end
        if ~isempty(RDS)
            save(stim_fname, '-append','RDS');
        end
        if ~isempty(RDS2)
            save(stim_fname, '-append','RDS2');
        end
        if ~isempty(RBS)
            save(stim_fname, '-append','RBS');
        end
        if ~isempty(RBS2)
            save(stim_fname, '-append','RBS2');
        end
        if ~isempty(NIS)
            save(stim_fname, '-append','NIS');
        end
        if ~isempty(RDSB)
            save(stim_fname, '-append','RDSB');
        end
        if ~isempty(RDSE)
            save(stim_fname, '-append','RDSE');
        end
    end


    function ExitStim
%         LoadIdentityClut(screen);
        Screen('CloseAll');
        Priority(0);
        ShowCursor;
    end

end


    
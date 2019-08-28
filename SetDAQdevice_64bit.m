% Set DAQ device for stim analog output and digital input from laser
% shutter.

% It is not possible to use the legacy interface on 64-bit Matlab
% (digitalio, putsample, etc. do not work!), but you must use the
% Session-Based Interface

if test == 0
    nidaq=[];
    if Param.useDAQdevice == 1
%       daq.getDevices  % to get info about installed devices
        daq.reset;
        nidaq = daq.createSession('ni');
        addAnalogOutputChannel(nidaq,'Dev1', 'ao0', 'Voltage');
        addDigitalChannel(nidaq,'Dev1', 'port0/line0', 'InputOnly');
        outputSingleScan(nidaq,0);
        inputSingleScan(nidaq);
    %     daqreset  % reset DAQ devices
    %     dio = digitalio('nidaq','Dev1');
    %     addline(dio,0,'in');     %Ale: input of the trigger (shutter); should be 0 if closed, 1 if open.
    %     ao=analogoutput('nidaq','Dev1');
    %     addchannel(ao,0);
    %     putsample(ao,0);
    elseif Param.useArduino == 1
        nidaq=[];
        Param.ArdPinNr = 9;
        global Ard %#ok<TLEV>
        if isempty(Ard)
            Ard = arduino('COM3');
%             Ard = arduino('/dev/ttyUSB0');
            Ard.pinMode(Param.ArdPinNr,'output');
        else
            disp('Arduino already connected !');
        end
        Ard.digitalWrite(Param.ArdPinNr,0); % set to 0V
    end
else
    nidaq=[];
%     ao=[];
%     dio=[];
end
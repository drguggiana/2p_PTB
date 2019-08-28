% set DAQ device for stim analog output

if test == 0 
    % these two global variables are used in callback function in udp client
    global mystring
    global stopRequest
    global stopScan
    mystring{1} = '';
    global trigger
    % global triggerRise
    daqreset  % reset DAQ devices
    dio = digitalio('nidaq','Dev1');
    addline(dio,0,'in');     %Ale: input of the trigger (shutter); should be 0 if closed, 1 if open.
    ao=analogoutput('nidaq','Dev1');
    addchannel(ao,0);
    putsample(ao,0);
    stopRequest = false;
else
    ao=[];
    dio=[];
end
function StimSignal( Param, test, nidaq, value )


if test == 0
    
    if     isfield(Param,'useDAQdevice')...
            && Param.useDAQdevice==1
        
        % putsample(ao,value); 
        outputSingleScan(nidaq,value);
        
    elseif isfield(Param,'useArduino')...
            && Param.useArduino==1
        
        global Ard %#ok<TLEV>
        if value > 0
            Ard.digitalWrite(Param.ArdPinNr,1); % set 5V
        else % value==0
            Ard.digitalWrite(Param.ArdPinNr,0); % set 0V
        end
        
    end
    
end
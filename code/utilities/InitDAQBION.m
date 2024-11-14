function dq = InitDAQBION
% create the daq object for reading and writing
% must not change unless hardware setup does

try
    dq = daq("ni");
    addinput(dq,"Dev1","port0/line0:7,port1/line0:7","Digital");
    addoutput(dq, "Dev1", "port3/line0:7,port4/line0:7", "Digital");
    write(dq, zeros(1,16));
catch
    error("this failed");
end


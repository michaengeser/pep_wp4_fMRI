function [TimeStamp, scantick] = GetTriggerDAQBION(daq_var, varargin)


scantick = 0;
while ~scantick
    D=read(daq_var,1,"OutputFormat","Matrix");
    if nargin > 1
        scantick = D(varargin{1});
    else
        scantick = D(16);
    end
end
TimeStamp = datetime;
%TimeStamp = GetSecs; % time stamp of the scan (Psychtoolbox)



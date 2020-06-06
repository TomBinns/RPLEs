function conv_times = convTimes(times,fs)
% convTimes - Converts times from mrk format to cnt format
% 
% Argument(s):
%  times -    times of markers from mrk (mrk.time)
%  fs -       sampling rate of the EEG data (cnt.fs)
% 
% Returns:
%  conv_times - times of markers in cnt format
% 
% Thomas Binns, 2020

si = 1000/fs;
timeeps = si/100;

conv_times = nan(size(times));
for aa = 1:length(times)
    conv_times(aa) = ceil((times(aa)-timeeps)/si);
end


end
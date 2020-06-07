function [mrk,idx] = removeShortTrials(mrk,minTime,classes)
% removeShortTrials - Removes trial markers for trials in which the 
%     time between two classes is shorter than a given time
% 
% Argument(s):
%  mrk -      trial markers
%  minTime -  minimum amount of acceptable time between markers of two 
%             classes in each trial
%  classes -  classes for which the time between is checked in each trial 
%             (e.g. trial start and movement onset)
% 
% Returns:
%  mrk -      new trial markers with short trials removed
%  idx -      indexes of removed trials, if any, else empty
% 
% Author(s): Thomas Binns, 2020


WT = getWTs(mrk,classes);

AcceptTrials = zeros(1,length(WT));
Acceptmrk = zeros(1,(length(mrk.time)));
bb = 0;
for aa = 1:length(WT)
    if WT(aa) > minTime
        AcceptTrials(aa)=1;
        for cc = 1:size(mrk.y,1)
            bb=bb+1;
            Acceptmrk(bb)=1;
        end
    else
        bb=bb+size(mrk.y,1);
    end
end

toDel = zeros(1,(length(Acceptmrk)-sum(Acceptmrk)));
dd = 1;
for cc = 1:length(Acceptmrk)
    if Acceptmrk(cc) == 0
        toDel(dd)=cc;
        dd = dd+1;
    end
end

mrk.y(:,toDel) = [];
mrk.time(toDel) = [];

idx = find(AcceptTrials==0);


end

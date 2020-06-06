function WT = getWTs(mrk,classes)
% getWTs - Gets waiting times between two classes in each trial
% 
% Argument(s):
%  mrk -      trial markers
%  classes -  classes for which the time between will be taken, for each 
%             trial (e.g. trial start and movement onset)
% 
% Returns:
%  WT -       times between the chosen classes for each trial
% 
% Author(s): Thomas Binns, 2019
 
mrk = mrk_selectClasses(mrk,classes);
convtimes = convTimes(mrk.time,100);
WT = convtimes(logical(mrk.y(2,:)))-convtimes(logical(mrk.y(1,:)));
WT = WT*10;

end
function centres = getErrorBarCentre(nGroups,nBars,group)
% Finds what the centres of error bars should be for stacked bar charts

% INPUTS:
% nGroups - number of groups that are being plotted
% nBars - number of bars being plotted in each group
% group - number of the current group being plotted

% OUTPUT:
% centres - locations for the centres of each error bar in this particular
%           group ([1 x nBars] row vector)


% find width for each bar group
groupWidth = min(0.8,nBars/(nBars+1.5));
centres = (1:nGroups) - groupWidth/2 + (2*group-1) * groupWidth / (2*nBars);

end
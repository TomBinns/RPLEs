function [rplemrk,rpleTs,NoRPLEs,lastadd] = findRPLEs(Cout,...
    threshold,mrk,untilTime,Cival)
% rple_findRPLEs - finds RPLEs in the EEG data
% 
% Argument(s):
%  Cout -      classifier output for the current subject for the current
%              mode
%  threshold - threshold for what counts as a RPLE
%  mrk -       trial markers
%  untilTime - time until which RPLEs are found in each trial; if 0, RPLEs
%              are found until movement onset
%  Cival -     Takes the form [X,Y]
%              X = time that must occur between RPLEs, threshold 
%              crossings that occur within this time of a previous 
%              threshold crossing are excluded, starting from the earliest 
%              crossing of the trial; if 0, no threshold crossings excluded
%              Y = end of the 'hit' window (the window where a threshold 
%              crossing counts as the classifier detecting a RP), if a 
%              crossing occurs within X of a crossing during the 'hit' 
%              window, this crossing outside the 'hit' window is removed; 
%              if NaN, this does not happen; 
%              'hit' window is specified as a negative value (i.e. -500 = 
%              500 ms before movement onset)
% 
% Returns:
%  rplemrk -   trial markers with RPLE markers added
%  rpleTs -    index of where RPLEs occured (which trials and at what 
%              times)
%  NoRPLEs -   if 0, the data contained RPLEs; if 1, the data contained no 
%              RPLEs
%  lastadd -   number of RPLEs in each trial
% 
% Author(s): Thomas Binns, 2019


if length(mrk.className) ~= 3
    error('This function is not designed for this mrk format');
end

% checks to see if any Couts are too short and do not contain the hit
% window; if so, excludes them
zz = 1;
shorttrials = zeros(1,length(Cout));
for aa = 1:length(Cout)
    if Cout{aa}.t(1) > untilTime
        shorttrials(zz) = aa;
        zz = zz+1;
    end
end
shorttrials(zz:end) = [];
if ~isempty(shorttrials)
    Cout(shorttrials) = [];
    zz = 1;
    for aa = 1:length(shorttrials)
        idx(zz:zz+2) = shorttrials(aa)*3-2:1:shorttrials(aa)*3;
        zz = zz+3;
    end
    mrk.y(:,idx) = [];
    mrk.time(idx) = [];
end
        
if Cival(1) ~= 0 && ~isnan(Cival(2))
    % find where the position of the hit window start is for each trial
    winendTimes = zeros(1,length(Cout));
    for aa = 1:length(Cout)
        winendTimes(aa) = find(Cout{aa}.t==Cival(2));
    end
    % find the longest trial so that the matrix for the times where Cout
    % changes can be sufficiently large
    Coutsizes = zeros(1,length(Cout));
    for aa = 1:length(Cout)
        Coutsizes(aa) = length(Cout{aa}.t);
    end
    CoutChangeTs = NaN(length(Cout),max(Coutsizes));
    % finds where Cout crosses the threshold
    for aa = 1:length(Cout)
        zz = 1;
        for bb = 1:winendTimes(aa)
            if Cout{aa}.x(bb) <= threshold && Cout{aa}.x(bb+1) >=...
                    threshold && Cout{aa}.t(bb+1) <= Cival(2)
                CoutChangeTs(aa,zz) = Cout{aa}.t(bb+1);
                zz = zz+1;
            end
        end
    end
    % cleans up CoutChangeTs so that no RPLEs occur within 1 classifier...
    % interval duration of each other removes all crossings that occur...
    % within 1 interval of the final crossing occuring during the hit...
    % period, starting with the last crossing and moving backwards, so...
    % that this is taken as the marker of the 'true RP'
    Ncrossings = zeros(length(Cout),1);
    for aa = 1:length(Cout)
        Ncrossings(aa) = (find(isnan(CoutChangeTs(aa,:)),1)-1);
        for bb = Ncrossings(aa):-1:1
            if bb - 1 > 0
                for cc = bb-1:-1:1
                    if CoutChangeTs(aa,bb) >= untilTime
                        if ~isnan(CoutChangeTs(aa,bb)) ||...
                                ~isnan(CoutChangeTs(aa,cc))
                            if CoutChangeTs(aa,bb) - CoutChangeTs(aa,cc)...
                                    <= Cival(1)
                                CoutChangeTs(aa,cc) = NaN;
                            end
                        end
                    end
                end
            end
        end
    end
    CoutChangeTs = sort(CoutChangeTs,2);
    % removes any crossings that occur within 1 interval of other...
    % crossings outside the hit period, starting with the first crossing...
    % and moving forwards
    for aa = 1:length(Cout)
        Ncrossings(aa) = (find(isnan(CoutChangeTs(aa,:)),1)-1);
        for bb = 1:Ncrossings(aa)
            for cc = bb+1:Ncrossings(aa)
                if ~isnan(CoutChangeTs(aa,bb)) ||...
                        ~isnan(CoutChangeTs(aa,cc))
                    if abs(CoutChangeTs(aa,bb) - CoutChangeTs(aa,cc))...
                            <= Cival(1)
                        CoutChangeTs(aa,cc) = NaN;
                    end
                end
            end
        end
    end
    CoutChangeTs = sort(CoutChangeTs,2);
    for aa = 1:length(Cout)
        Ncrossings(aa) = (find(isnan(CoutChangeTs(aa,:)),1)-1);
        for bb = 1:Ncrossings(aa)
            if CoutChangeTs(aa,bb) >= untilTime
                CoutChangeTs(aa,bb) = NaN;
            end
        end
    end
    CoutChangeTs = sort(CoutChangeTs,2);
end

if isnan(Cival(2))
    % find where the position of the hit window start is for each trial
    winstartTimes = zeros(1,length(Cout));
    for aa = 1:length(Cout)
        winstartTimes(aa) = (find(Cout{aa}.t==untilTime));
    end
    % find the longest trial so that the matrix for the times where Cout...
    % changes can be sufficiently large for all possibilities
    Coutsizes = zeros(1,length(Cout));
    for aa = 1:length(Cout)
        Coutsizes(aa) = length(Cout{aa}.t);
    end
    CoutChangeTs = NaN(length(Cout),max(Coutsizes));
    % finds where Cout crosses the threshold
    for aa = 1:length(Cout)
        zz = 1;
        for bb = 1:winstartTimes(aa)
            if Cout{aa}.x(bb) <= threshold &&...
                    Cout{aa}.x(bb+1) >= threshold &&...
                    Cout{aa}.t(bb+1) <= untilTime
                CoutChangeTs(aa,zz) = Cout{aa}.t(bb+1);
                zz = zz+1;
            end
        end
    end
end
    
if Cival(1) ~= 0 && isnan(Cival(2))
    % removes any crossings that occur within 1 interval of other...
    % crossings outside the hit period, starting with the first...
    % crossing and moving forwards
    Ncrossings = zeros(1,length(Cout));
    for aa = 1:length(Cout)
        Ncrossings(aa) = find(isnan(CoutChangeTs(aa,:)),1)-1;
        for bb = 1:Ncrossings(aa)
            for cc = bb+1:Ncrossings(aa)
                if ~isnan(CoutChangeTs(aa,bb)) ||...
                        ~isnan(CoutChangeTs(aa,cc))
                    if abs(CoutChangeTs(aa,bb) - CoutChangeTs(aa,cc))...
                            <= Cival(1)
                        CoutChangeTs(aa,cc) = NaN;
                    end
                end
            end
        end
    end
    CoutChangeTs = sort(CoutChangeTs,2);
end

% Find the times (in the form used by mrk.time) at which crossings occur...
% and how many times this occurs per trial
yy = 0;
zz = 1;
newtimes = zeros(1,max(Coutsizes)*length(Cout));
lastadd = zeros(1,length(Cout));
for aa = 1:length(Cout)
    lastadd(aa) = find(isnan(CoutChangeTs(aa,:)),1)-1;
    for bb = 1:lastadd(aa)
        if ~isnan(CoutChangeTs(aa,bb))
            newtimes(zz) = mrk.time(yy+2)+CoutChangeTs(aa,bb);
            zz = zz+1;
        end
    end
    yy = yy+3;
end
newtimes(newtimes==0) = [];

if sum(lastadd) == 0
    rplemrk = mrk;
    NoRPLEs = 1;
else
    NoRPLEs = 0;
    % create a new mrk that will contain the RPLE markers
    rplemrk = mrk;
    % create new class names for the RPLE markers
    rplemrk.className = {'trial start','rple','movement onset',...
        'trial end'};
    % increase the size of the marker matrix to account for the additional...
        % RPLE markers
    rplemrk.y = zeros(length(rplemrk.className),...
        (length(mrk.time)+sum(lastadd)));
    % find trials in which there are no RPLEs
    norples = find(lastadd(1,:)==0);
    % add the new times where RPLEs occur and arrange them into the...
    % correct order based on time
    rplemrk.time = cat(2,rplemrk.time,newtimes);
    rplemrk.time = sort(rplemrk.time);
    % fills in the marker (.y) field with the correct markers in the...
    % correct places
    zz = 1;
    cc = 1;
    for aa = 1:length(Cout)
        % adds the trial start marker
        rplemrk.y(1,cc) = mrk.y(1,zz);
        zz = zz+1;
        cc = cc+1;
        % if the trial contained at least one RPLE, add the appropriate...
        % number of RPLE markers to the RPLE rows
        if lastadd(aa) > 0 % if the trial contain at least one RPLE, add...
            % the appropriate no. of RPLE markers to the RPLE rows
            rplemrk.y(2,cc:cc+lastadd(aa)-1) = 1;
            cc = cc+lastadd(aa);
        else % otherwise, leave the RPLE row blank and move on to the next row
            cc = cc+1;
        end
        % adds the movement onset marker
        rplemrk.y(end-1,cc) = mrk.y(2,zz);
        zz = zz+1;
        cc = cc+1;
        % adds the trial end marker
        rplemrk.y(end,cc) = mrk.y(3,zz);
        zz = zz+1;
        cc = cc+1;
    end
    % adds time points to the .t field (for trials in which there are no...
    % RPLEs) so that the time points can be arranged based on value and...
    % the markers will still correspond to the correct time points
    timesToAdd = zeros(1,length(norples));
    for aa = 1:length(norples)
        % for trials in which there are no RPLEs, add markers which will...
        % occur after trial start
        timesToAdd(aa) = mrk.time((norples(aa)*3)-2)+(.1);
    end
    rplemrk.time = cat(2,rplemrk.time,timesToAdd);
    rplemrk.time = sort(rplemrk.time);
end

rpleTs = repmat({struct('rpleTs',{})},size(CoutChangeTs,1),1);
for aa = 1:size(CoutChangeTs,1)
    rpleTs{aa} = CoutChangeTs...
        (aa,1:find(isnan(CoutChangeTs(aa,:)),1)-1);
end

end

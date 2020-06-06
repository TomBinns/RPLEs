function [thresholds,SecsPerEvent] = findThresholds(Cout_p1)
% findThresholds - Generates thresholds 1 - 5 to apply to the EEG data to 
%     identify RPLEs
%     
% Argument(s):
%  Cout_p1 -  classifier outputs of phase 1 EEG data of all subjects and 
%             each mode (spatio-temporal, temporal, and spatial)
% 
% Returns:
%  thresholds - thresholds 1 - 5 for each subject and all three modes
%  SecsPerEvent - number of seconds per event for threshold 5, for all 
%                 subjects and all modes (should be close to 8.95)
% 
% Author(s): Thomas Binns, 2020                


global opt

thresholds = cell(1,3);
for aa = 1:3
    thresholds{aa} = nan(5,length(opt.subjs_all));
end

T5acc = nan(length(opt.subjs_all),3);
for aa = 1:3
    if aa == 1
        mode = 'ST';
    elseif aa == 2
        mode = 'T';
    else
        mode = 'S';
    end
    P1Cout = Cout_p1(:,aa);
    for bb = 1:length(opt.subjs_all)
        subj_code = opt.subjs_all{bb};
        % loads phase 1 data
        [mrk,cnt] = loadData(subj_code,'Phase1');
        % selects trials with movement onset
        mrk = mrk_selectClasses(mrk,{'trial start','movement onset',...
            'trial end'});
        trial_mrk = getTrialMarkers(mrk);
        trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
        mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
        % removes trials in which no RPLEs can occur due to the time lost
        % resulting from the classifier training and the RP window
        [mrk,idx] = removeShortTrials(mrk,(1000+abs(opt.untilTime)),...
            {'trial start','movement onset'});
        P1Cout{bb}(idx) = [];
        % gets waiting times of trials
        WT = getWTs(mrk,{'trial start','movement onset'});
        WT = WT + opt.untilTime - 1000;
        WT = WT/1000;
        meanWT = mean(WT);
        % gets threshold for classifier output
        edges = [-Inf opt.WindowStartPoint opt.WindowEndPoint Inf];
        opt.method = 'setval';
        [highthresh,NRPLEs] = getThreshMode(mode,P1Cout{bb},edges,...
            mrk,cnt,meanWT);
        T5acc(bb,aa) = (NRPLEs/sum(WT))*8.95;
        opt.method = 'avgcout';
        lowthresh = getThreshMode(mode,P1Cout{bb},edges,mrk,cnt,meanWT);
        thresholds{aa}(:,bb) = linspace(lowthresh,highthresh,5);
    end
end

SecsPerEvent = nan(size(T5acc));
SecsPerEvent(:) = 8.95;
SecsPerEvent = SecsPerEvent./T5acc;

end
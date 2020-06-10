function [rpleepos,Eps,rpleTs] = modeAnalysisCz(P1Cout,...
    RTCout,varargin)
% modeAnalysisCz - Finds RPLEs in EEG data for the Cz-only mode
% 
% Argument(s):
%  P1Cout -   classifier output (Cz-only) for phase 1
%  RTCout -   classifier output (Cz-only) for phase RT
%  OPTIONAL ARGUMENT(S):
%   'thresh' - name-value pair; if present, the following input is taken as
%              the thresholds
%   'avg' -    name-value pair; specifies whether the epoched data should
%              be averaged (1, default) or not (0)
%   'selMode' - name-value pair; if present, the following input specifies
%               which modes the analysis should be conducted on; takes the 
%               form [spatio-temporal,temporal,spatial], with 1 signalling 
%               the mode should be analysed, and 0 signalling the mode 
%               should be skipped (e.g. [1,0,1] = spatio-temporal and 
%               spatial mode analysed, temporal skipped)
% 
% Returns:
%  rpleepos - epoched data containing RPs and RPLEs
%  Eps -      number of events per second of analysed data for each mode,
%             threshold, and subject
%  rpleTs -   indexes of when RPLEs occured for each mode, threshold, and
%             subject
%  rsqs -     r^2 values for the spatial mode of analysis scalp plots
%  thresholds - thresholds used to find RPLEs
% 
% Author(s): Thomas Binns, 2020

global opt
warning off

avg = 1;
% load optional arguments
if ~isempty(varargin)
    for aa = 1:length(varargin)/2
        switch varargin{(aa*2)-1}
            case 'thresh'
                allthresh = varargin{aa*2};
            case 'avg'
                avg = varargin{aa*2};
            otherwise
                error('Unexpected input argument');
        end
    end
end

rpleepos = repmat({struct('rpleepos',{})},1,1);
rpleTs = repmat({struct('rpleTs',{})},1,1);
Eps = repmat({struct('Eps',{})},1,1);
for aa = 1:1
    rpleepos{aa} = repmat({struct('rpleepos',{})},1,length(opt.subjs_all));
    rpleTs{aa} = repmat({struct('rpleTs',{})},2,length(opt.subjs_all));
    Eps{aa} = nan(2,length(opt.subjs_all));
end


P1Cout1 = P1Cout(:,1);
RTCout1 = RTCout(:,1);

for aa = 1:length(opt.subjs_all)
    subj_code = opt.subjs_all{aa};
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
    P1Cout1{aa}(idx) = [];
    % gets waiting times of trials
    WT = getWTs(mrk,{'trial start','movement onset'});
    WT = WT + opt.untilTime - 1000;
    WT = sum(WT)/1000;
    % gets threshold for classifier output
    threshold = allthresh{1}(aa);
    % Gets info about RPLEs and the times of threshold crossings
    [rplemrk,CoutChangeTs,NoRPLEs,lastadd] = findRPLEs(P1Cout1{aa},...
        threshold,mrk,opt.untilTime,opt.Cival);
    rpleTs{1}{1,aa} = CoutChangeTs;
    % Isolates RPs & RPLEs and gets number of RPLEs per second
    if NoRPLEs == 0 % if RPLEs present
        rplemrk = mrk_selectClasses(rplemrk,{'rple','movement onset'});
        rplemrk2 = mrk_selectClasses(rplemrk,{'rple'});
        epo = proc_segmentation(cnt,rplemrk2,opt.epochSegment);
        epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
        [~,iArte] = proc_rejectArtifactsMaxMin(epo,200);
        Eps{1}(1,aa) = (sum(lastadd)-length(iArte))/WT;
    else % if no RPLEs, just epoch RPs
        rplemrk = mrk_selectClasses(rplemrk,{'movement onset'});
        Eps{1}(1,aa) = 0;
    end
    % Epochs data and performs artefact rejection
    epo = proc_segmentation(cnt,rplemrk,opt.epochSegment);
    epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
    epo.className = {'rple p1','movement onset'};
    markers = epo.y;
    [epo,iArte] = proc_rejectArtifactsMaxMin(epo,200);
    rpleepos{1}{1,aa} = epo;
    % removes RPLE markers (from rpleTs) for events that are rejected
    iArte(markers(2,iArte)==1) = [];
    for bb = 1:length(iArte)
        iArte(bb) = iArte(bb) - sum(markers(2,1:iArte(bb)));
    end
    zz = 1;
    for bb = 1:length(rpleTs{1}{aa})
        if ~isempty(rpleTs{1}{aa}{bb})
            for cc = 1:length(rpleTs{1}{aa}{bb})
                if ismember(zz,iArte)
                    rpleTs{1}{aa}{bb}(cc) = nan;
                end
                zz = zz+1;
            end
            rpleTs{1}{aa}{bb}(isnan(rpleTs{1}{aa}{bb})) = [];
        end
    end
    
    % loads phase RT data
    [mrk,cnt] = loadData(subj_code,'RT');
    % selects trials with movement onset
    mrk = mrk_selectClasses(mrk,{'trial start','go signal','trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    % removes trials in which no RPLEs can occur due to the time lost
    % resulting from the classifier training and the RP window
    [mrk,idx] = removeShortTrials(mrk,1000,...
        {'trial start','go signal'});
    RTCout1{aa}(idx) = [];
    % gets waiting times of trials
    WT = getWTs(mrk,{'trial start','go signal'});
    WT = WT - 1000;
    WT = sum(WT)/1000;
    % Gets info about RPLEs and the times of threshold crossings
    [rplemrk,CoutChangeTs,NoRPLEs,lastadd] = findRPLEs(RTCout1{aa},...
        threshold,mrk,0,opt.Cival);
    rpleTs{1}{2,aa} = CoutChangeTs;
    % Isolates RPs & RPLEs and gets number of RPLEs per second
    if NoRPLEs == 0 % if RPLEs present
        rplemrk = mrk_selectClasses(rplemrk,{'rple'});
        epo = proc_segmentation(cnt,rplemrk,opt.epochSegment);
        epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
        [~,iArte] = proc_rejectArtifactsMaxMin(epo,200);
        Eps{1}(2,aa) = (sum(lastadd)-length(iArte))/WT;
        % Epochs data and performs artefact rejection
        epo = proc_segmentation(cnt,rplemrk,opt.epochSegment);
        epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
        epo.className = {'rple rt'};
        [epo,iArte] = proc_rejectArtifactsMaxMin(epo,200);
        % removes RPLE markers (from rpleTs) for events that are rejected
        zz = 1;
        for bb = 1:length(rpleTs{1}{aa})
            if ~isempty(rpleTs{1}{aa}{bb})
                for cc = 1:length(rpleTs{1}{aa}{bb})
                    if ismember(zz,iArte)
                        rpleTs{1}{aa}{bb}(cc) = nan;
                    end
                    zz = zz+1;
                end
                rpleTs{1}{aa}{bb}(isnan(rpleTs{1}{aa}{bb})) = [];
            end
        end
        % combines phase RT RPLEs into the same variable as the phase 1 RPs
        % and RPLEs
        rpleepos{1}{1,aa}.x = cat(3,rpleepos{1}{1,aa}.x,epo.x);
        rpleepos{1}{1,aa}.y(3,:) = 0;
        rty = zeros(3,size(epo.x,3));
        rty(3,:) = 1;
        rpleepos{1}{1,aa}.y = cat(2,rpleepos{1}{1,aa}.y,rty);
        rpleepos{1}{1,aa}.className = cat(2,rpleepos{1}{1,aa}.className,...
            epo.className);
        if avg == 1
            rpleepos{1}{1,aa} = proc_average(rpleepos{1}{1,aa},'Stats',1);
        end
    else % if no RPLEs, set number of events to 0 and use previous epochs
             % as a placeholder
        Eps{1}(2,aa) = 0;
        if avg == 1
            rpleepos{1}{1,aa} = proc_average(rpleepos{1}{1,aa});
        end
    end
end


end
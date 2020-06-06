function [threshold,numRPLEs] = getThreshMode(mode,Cout,edges,mrk,...
    cnt,meanWT)
% getThreshMode - Gets the threshold for the classifier output to
%     identify RPLEs
% 
% Argument(s):
%  mode -     mode for which the threshold is being generated; ST
%             (spatio-temporal), T (temporal), or S (spatial)
%  Cout -     classifier output of the subject for the current mode
%  edges -    window during which a threshold crossing would count as the
%             classifier correctly detecting a RP
%  mrk -      trial markers
%  cnt -      EEG data
%  meanWT -   mean waiting time from trial start to movement onset for the
%             subject
% 
% Returns:
%  threshold - threshold for which a crossing indicates the occurence of a
%              RPLE
%  numRPLEs - number of RPLEs identified by this threshold
% 
% Author(s): Thomas Binns, 2020

global opt

if ~strcmp(opt.method,'setval') && ~strcmp(opt.method,'avgcout') &&...
        ~strcmp(opt.method,'default')
    error('The requested method is not recognised');
end

if ~strcmp(mode,'ST') && ~strcmp(mode,'T') && ~strcmp(mode,'S')
    error('The requested mode is not recognised');
end

numRPLEs = [];
if strcmp(opt.method,'setval')
    % gets threshold based on whichever threshold gives the closest number
    % of crossings to the target
    
    % gets the regular threshold using the F-score in case this is needed
    % as a tie-breaker
    x_all = cellfun(@(f)f.x,Cout,'UniformOutput',false);
    for aa = 1:length(x_all)
        x_all{aa} = transpose(x_all{aa});
    end
    x_all = [x_all{:}];
    thresh = linspace(0,prctile(x_all,99.5),100);
    Nth = length(thresh);
    Nt = length(Cout);
    T = inf(Nt,Nth);
    for jj = 1:Nt
        tind = Cout{jj}.t<=edges(3);
        x = Cout{jj}.x(tind);
        t = Cout{jj}.t(tind);
        for kk = 1:Nth
            ind = find(diff(sign(x-thresh(kk)))==2,1);
            if not(isempty(ind))
                T(jj,kk) = t(ind);
            end
        end
    end
    R = zeros(Nth,length(edges)-1);
    for kk = 1:Nth
        R(kk,:) = histcounts(T(:,kk),edges,'normalization','probability');
    end
    F = fScore(R(:,1),R(:,2),R(:,3),opt.beta);
    F = smooth(F,10);
    [~,ind_maxF] = max(F);
    
    % looks at the number of crossings (after artefact exclusion) that each
    % threshold gives
    NRPLEs = nan(1,length(thresh));
    for aa = 1:length(thresh)
%         [rplemrk,~,~,lastadd] = rple_findRPLEs(Cout,thresh(aa),mrk,...
%             opt.untilTime,opt.Cival,'cat');
        [rplemrk,~,~,lastadd] = findRPLEs(Cout,thresh(aa),mrk,...
            opt.untilTime,opt.Cival);
        NRPLEs(aa) = sum(lastadd);
        if NRPLEs(aa) > 0
            rplemrk = mrk_selectClasses(rplemrk,'rple');
            if strcmp(mode,'ST')
                epo = proc_segmentation(cnt,rplemrk,opt.epochSegment);
                epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
            elseif strcmp(mode,'T')
                epo = proc_segmentation(cnt,rplemrk,opt.epochSegment);
                epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
                epo = proc_meanAcrossChannels(epo);
            elseif strcmp(mode,'S')
                epo = proc_segmentation(cnt,rplemrk,[-400,0]);
                epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
                epo = proc_meanAcrossTime(epo);
            end
            [~,iArte] = proc_rejectArtifactsMaxMin(epo,200);
            NRPLEs(aa) = NRPLEs(aa)-length(iArte);
        end
    end
    
    % find which threshold gives the number of RPLEs that is closest to the
    % target
    perSecTarget = 8.95;
    perTrialTarget = meanWT/perSecTarget;
    target = length(Cout)*perTrialTarget;
%     target = length(Cout)/2;
    ithresh = find(NRPLEs==target);
    if ~isempty(ithresh)
        if length(ithresh) > 1
            [~,ind] = min(abs(ind_maxF-ithresh'));
            ithresh = ithresh(ind);
        end        
    else
        [~,ithresh] = min(abs(target-NRPLEs'));   
    end
    threshold = thresh(ithresh);
    numRPLEs = NRPLEs(ithresh);
    
elseif strcmp(opt.method,'avgcout')
    % Gets threshold based on average classifier output at time of movement
    % onset
    Coutvals = nan(1,length(Cout));
    for aa = 1:length(Cout)
        iCout = find(Cout{aa}.t==0);
        Coutvals(aa) = Cout{aa}.x(iCout);
    end
    threshold = mean(Coutvals);
else
     x_all = cellfun(@(f)f.x,Cout,'UniformOutput',false);
    for aa = 1:length(x_all)
        x_all{aa} = transpose(x_all{aa});
    end
    x_all = [x_all{:}];
    thresh = linspace(0,prctile(x_all,99.5),100);
    Nth = length(thresh);
    Nt = length(Cout);
    T = inf(Nt,Nth);
    for jj = 1:Nt
        tind = Cout{jj}.t<=edges(3);
        x = Cout{jj}.x(tind);
        t = Cout{jj}.t(tind);
        for kk = 1:Nth
            ind = find(diff(sign(x-thresh(kk)))==2,1);
            if not(isempty(ind))
                T(jj,kk) = t(ind);
            end
        end
    end
    R = zeros(Nth,length(edges)-1);
    for kk = 1:Nth
        R(kk,:) = histcounts(T(:,kk),edges,'normalization','probability');
    end
    F = fScore(R(:,1),R(:,2),R(:,3),opt.beta);
    F = smooth(F,10);
    [~,ind_maxF] = max(F);
    
    threshold = thresh(ind_maxF);
end

end
function [rpleepos,Eps,rpleTs] = modeAnalysisNoise...
    (Cout_noise,thresholds,type)
% modeAnalysisNoise - Finds RPLEs in EEG data for the Cz-only mode
% 
% Argument(s):
%  Cout_noise - classifier output (Cz-only) for noise data
%  thresholds - thresholds to find RPLEs
%  type -     noise data which should be loaded; string of the form
%             'NoiseXXX' where XXX specify the alpha value of the noise
%             (e.g. 000 = 0.00; 150 = 1.50; Per = Personal)
% 
% Returns:
%  rpleepos - epoched data containing RPs and RPLEs
%  Eps -      number of events per second of analysed data for each mode,
%             threshold, and subject
%  rpleTs -   indexes of when RPLEs occured for each mode, threshold, and
%             subject
% 
% Author(s): Thomas Binns, 2020


global opt
warning on

rpleepos = repmat({struct('rpleepos',{})},1,1);
rpleTs = repmat({struct('rpleTs',{})},1,1);
Eps = repmat({struct('Eps',{})},1,1);
for aa = 1:1
    rpleepos{aa} = repmat({struct('rpleepos',{})},1,length(opt.subjs_all));
    rpleTs{aa} = repmat({struct('rpleTs',{})},1,length(opt.subjs_all));
    Eps{aa} = nan(1,length(opt.subjs_all));
end


for aa = 1:length(opt.subjs_all)
    subj_code = opt.subjs_all{aa};
    [mrk,cnt] = loadData(subj_code,type);
    % gets WT of data
    WT = getWTs(mrk,{'trial start','movement onset'});
    WT = WT - 1000;
    WT = sum(WT)/1000;
    % Gets info about RPLEs and the times of threshold crossings
    [rplemrk,CoutChangeTs,NoRPLEs,lastadd] = findRPLEs(Cout_noise{aa},...
        thresholds(aa),mrk,0,opt.Cival);
    % Isolates RPs & RPLEs and gets number of RPLEs per second
    if NoRPLEs == 0 % if RPLEs present
        rplemrk = mrk_selectClasses(rplemrk,{'rple'});
        epo = proc_segmentation(cnt,rplemrk,opt.epochSegment);
        epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
        [~,iArte] = proc_rejectArtifactsMaxMin(epo,200);
        Eps{1}(aa) = (sum(lastadd)-length(iArte))/WT; % number of events
            % per second
        % Repackages CoutChangeTs into cells and removes excessive NaNs
        rpleTs{1}{aa} = repmat({struct('rpleTs',{})},size(CoutChangeTs,1),1);
        for bb = 1:size(CoutChangeTs,1)
            rpleTs{1}{aa}{bb} = CoutChangeTs...
                (bb,1:find(isnan(CoutChangeTs(bb,:)),1)-1);
        end
        % Epochs data and performs artefact rejection
        epo = proc_segmentation(cnt,rplemrk,opt.epochSegment);
        epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
        epo.className = {'rple noise'};
        [epo,iArte] = proc_rejectArtifactsMaxMin(epo,200);
        % removes artefact-rejected RPLEs from rpleTs
        if ~isempty(iArte)
            zz = 1;
            for bb = 1:size(rpleTs{1}{aa},1)
                if isempty(rpleTs{1}{aa}{bb})
                    zz = zz+1;
                else
                    for cc = 1:length(rpleTs{1}{aa}{bb})
                        if ismember(zz,iArte)
                            rpleTs{1}{aa}{bb}(zz) = nan;
                            zz = zz+1;
                        end
                    end
                end
                
            end
            for bb = 1:size(rpleTs{1}{aa},1)
                if ~isempty(rpleTs{1}{aa}{bb})
                    rpleTs{1}{aa}{bb}(isnan(rpleTs{1}{aa}{bb})) = [];
                end
            end
        end
%         rpleepos{1}{1,aa} = proc_average(epo,'Stats',1);
        rpleepos{1}{1,aa} = epo;
    else
        % if no RPLEs in the data, sets the number of events to 0 and
        % generates a placeholder epoch data struct
        Eps{1}(aa) = 0;
%         rpleepos{1}{1,aa} = proc_average(epo,'Stats',1);
        
        [mrk,cnt] = loadData(subj_code,'Phase1');
        mrk = mrk_selectClasses(mrk,{'trial start','movement onset',...
            'trial end'});
        trial_mrk = getTrialMarkers(mrk);
        trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
        mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
        mrk = mrk_selectClasses(mrk,'movement onset');
        epo = proc_segmentation(cnt,mrk,opt.epochSegment);
        epo = proc_baseline(epo,opt.baseln_len,opt.baseln_pos);
        
        epo.x = epo.x(:,:,1);
        epo.x(:,:) = nan;
        epo.y = 1;
        rpleepos{1}{1,aa} = epo;
    end
end


end
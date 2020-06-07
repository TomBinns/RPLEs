function [allDists,tables] = plotNoiseAnalysis(sets,thresholds)
% plotNoiseAnalysis - Plots results for the noise data (number of RPLEs,
%     similarity scores of RPLEs, appearance of RPLEs)
%     
% Argument(s):
%  sets -     contains the epoched RPLEs as well as information such as the 
%             number of RPLEs in the data for all subjects, for each 
%             threshold
%  thresholds - thresholds used to identify RPLEs
%  
% Returns:
%  allDists - similarity scores of the RPLEs for all subjects, and all 
%             thresholds
%  tables -   results tables for the statistical tests:
%              ttest2 - compares number of RPLEs in the EEG data to the
%                       noise data
%              kstest2 and ranksum - compare the similarity scores of EEG 
%                                    RPLEs to the noise RPLEs
%  
% Author(s): Thomas Binns, 2020 


global opt

%% plots number of events per second for noise data

Eps = cell(1,5);
for cc = 1:7
    if cc == 1
        for aa = 1:5
            Eps{aa} = nan(7,length(opt.subjs_all));
            for bb = 1:length(opt.subjs_all)
                p1Pos = find(strcmp(sets{aa}.rpleepos{1}{bb}.className,...
                    'rple p1'));
                rtPos = find(strcmp(sets{aa}.rpleepos{1}{bb}.className,...
                    'rple rt'));
                combWT = (sets{aa}.rpleepos{1}{bb}.N(p1Pos)/...
                    sets{aa}.Eps{1}(1,bb)) + ...
                    (sets{aa}.rpleepos{1}{bb}.N(rtPos)/...
                    sets{aa}.Eps{1}(2,bb));
                Eps{aa}(cc,bb) = (sets{aa}.rpleepos{1}{bb}.N(p1Pos)+...
                    sets{aa}.rpleepos{1}{bb}.N(rtPos))/combWT;
            end
        end
    else
        for aa = 1:5
            Eps{aa}(cc,:) = sets{aa}.noiseEps{cc-1};
        end
    end
end

SEM = cell(1,5);
for bb = 1:5
    for cc = 1:7
        SEM{bb}(1,cc) = std(Eps{bb}(cc,:))/sqrt(size(Eps{bb},2));
    end
end

figure
hold on
Colours = [.5,.5,.5;0,.45,.74;.93,.69,.13;.48,.18,.56;.3,.75,.93;...
    .47,.67,.19;.64,.08,.18;.85,.33,.1];
h = bar([1,2,3,4,5],[mean(Eps{1},2),mean(Eps{2},2),mean(Eps{3},2),...
    mean(Eps{4},2),mean(Eps{5},2)],'LineWidth',1.5);
h(1).Parent.FontSize = 20;
h(1).Parent.FontWeight = 'bold';
for aa = 1:7
    h(aa).FaceColor = Colours(aa,:);
    h(aa).Parent.XColor = [0,0,0];
    h(aa).Parent.YColor = [0,0,0];
end
for bb = 1:7
    centre = getErrorBarCentre(5,7,bb);
    errorbar(centre,[mean(Eps{1}(bb,:)),mean(Eps{2}(bb,:)),...
        mean(Eps{3}(bb,:)),mean(Eps{4}(bb,:)),mean(Eps{5}(bb,:))],...
        [SEM{1}(bb),SEM{2}(bb),SEM{3}(bb),SEM{4}(bb),SEM{5}(bb)],...
        'k','linestyle','none','LineWidth',3,'CapSize',10,...
        'Color',[0,0,0]);
end
xticks([1,2,3,4,5]);
xticklabels({'1','2','3','4','5'});
xlabel('Threshold')
ylabel('RPLEs/s');
legEntries = {'EEG','\alpha = 0','\alpha = 0.5','\alpha = 1',...
    '\alpha = 1.5','Personal \alpha','\alpha = 1.75'};
legend(legEntries,'Location','northeast');

% checks for sig diff
tables.tables = cell(1,3);
varNames = {'alpha = 0','alpha = 0.5','alpha = 1.0',...
    'alpha = 1.5','Personal alpha','alpha = 1.75'};
sigdiff = nan(5,6);
for bb = 1:6
    for aa = 1:5
        [~,sigdiff(aa,bb)] = ttest2(Eps{aa}(1,:),Eps{aa}(bb+1,:));
    end
end
tables.tables{1} = array2table(sigdiff,'VariableNames',varNames);
tables.desc{1} = 'ttest2: is there a significant difference between the number of RPLEs/s in the EEG data and the noise data?';


%% Gets combined EEG RPLE data

CoutDir = fullfile(opt.dataDir,'Cout');
loadname = fullfile(CoutDir,'Cout_Cz_p1.mat');
load(loadname);
loadname = fullfile(CoutDir,'Cout_Cz_rt.mat');
load(loadname);

apps = {'000','050','100','150','Per','175'};
dists = cell(1,5);
allDists = dists;
for aa = 1:5
    dists{aa} = cell(1,length(opt.subjs_all));
    allDists{aa}.rp = [];
    allDists{aa}.rple_eeg = [];
    allDists{aa}.rple_noise = cell(1,7);
    selthresh{1} = thresholds(aa,:);
    rpleepos = modeAnalysisCz(Cout_p1,Cout_rt,'thresh',selthresh,'noAvg',1);
    for bb = 1:length(rpleepos{1})
        rpleepos{1}{bb}.y(p1Pos,...
            logical(rpleepos{1}{bb}.y(rtPos,:))) = 1;
        rpleepos{1}{bb}.y(rtPos,:) = [];
        rpleepos{1}{bb}.className = {'rple eeg','movement onset'};
        sets{aa}.combepos{1}{bb} = rpleepos{1}{bb};
        sets{aa}.combepos{1}{bb}.y = cat(1,sets{aa}.combepos{1}{bb}.y,...
            zeros(6,size(sets{aa}.combepos{1}{bb}.y,2)));
        for cc = 1:6
            sets{aa}.combepos{1}{bb}.x = cat(3,sets{aa}.combepos{1}{bb}.x,...
                sets{aa}.noiseepos{cc}{bb}.x);
            newY = zeros(8,size(sets{aa}.noiseepos{cc}{bb}.y,2));
            newY(cc+2,:) = 1;
            sets{aa}.combepos{1}{bb}.y = cat(2,sets{aa}.combepos{1}{bb}.y,...
                newY);
            sets{aa}.combepos{1}{bb}.className{cc+2} = ...
                sprintf('rple noise %s',apps{cc});
        end
        % gets RP template
        eegPos = find(strcmp(sets{aa}.combepos{1}{bb}.className,...
            'rple eeg'));
        rpPos = find(strcmp(sets{aa}.combepos{1}{bb}.className,...
            'movement onset'));
        avgepo = proc_selectClasses(sets{aa}.combepos{1}{bb},...
            'movement onset');
        avgepo = proc_average(avgepo);
        template = avgepo.x;
        dists{aa}{bb}.rp = ...
            nan(sum(sets{aa}.combepos{1}{bb}.y(rpPos,:)),1);
        dists{aa}{bb}.rple_eeg = ...
            nan(sum(sets{aa}.combepos{1}{bb}.y(eegPos,:)),1);
        dists{aa}{bb}.rple_noise = cell(1,6);
        for cc = 1:6
            dists{aa}{bb}.rple_noise{cc} = ...
                nan(sum(sets{aa}.combepos{1}{bb}.y(cc+2,:)),1);
        end
        
        % applies RP template to single-trial RPs and RPLEs (EEG and noise)
        xx = 1;
        yy = 1;
        zz = ones(1,6);
        for cc = 1:size(sets{aa}.combepos{1}{bb}.y,2)
            event = sets{aa}.combepos{1}{bb}.x(:,:,cc);
            if ~isnan(norm(template-event))
                if sets{aa}.combepos{1}{bb}.y(rpPos,cc) == 1
                    dists{aa}{bb}.rp(xx) = norm(template-event);
                    xx = xx+1;
                elseif sets{aa}.combepos{1}{bb}.y(eegPos,cc) == 1
                    dists{aa}{bb}.rple_eeg(yy) = norm(template-event);
                    yy = yy+1;
                elseif sets{aa}.combepos{1}{bb}.y(3,cc) == 1
                    dists{aa}{bb}.rple_noise{1}(zz(1)) = norm(template-event);
                    zz(1) = zz(1)+1;
                elseif sets{aa}.combepos{1}{bb}.y(4,cc) == 1
                    dists{aa}{bb}.rple_noise{2}(zz(2)) = norm(template-event);
                    zz(2) = zz(2)+1;
                elseif sets{aa}.combepos{1}{bb}.y(5,cc) == 1
                    dists{aa}{bb}.rple_noise{3}(zz(3)) = norm(template-event);
                    zz(3) = zz(3)+1;
                elseif sets{aa}.combepos{1}{bb}.y(6,cc) == 1
                    dists{aa}{bb}.rple_noise{4}(zz(4)) = norm(template-event);
                    zz(4) = zz(4)+1;
                elseif sets{aa}.combepos{1}{bb}.y(7,cc) == 1
                    dists{aa}{bb}.rple_noise{5}(zz(5)) = norm(template-event);
                    zz(5) = zz(5)+1;
                elseif sets{aa}.combepos{1}{bb}.y(8,cc) == 1
                    dists{aa}{bb}.rple_noise{6}(zz(6)) = norm(template-event);
                    zz(6) = zz(6)+1;
                else
                    error('An event marker is missing')
                end
            end
        end
        allDists{aa}.rp = cat(1,allDists{aa}.rp,dists{aa}{bb}.rp);
        allDists{aa}.rple_eeg = cat(1,allDists{aa}.rple_eeg,...
            dists{aa}{bb}.rple_eeg);
        for cc = 1:6
            allDists{aa}.rple_noise{cc} = cat(1,allDists{aa}.rple_noise{cc},...
                dists{aa}{bb}.rple_noise{cc});
        end
    end
end

for aa = 1:5
    for bb = 1:6
        if isnan(sum(allDists{aa}.rple_noise{bb}))
            idx = find(isnan(allDists{aa}.rple_noise{bb}));
            allDists{aa}.rple_noise{bb}(idx) = [];
        end
    end
end

% gets and plots medians for data of individual participants
indivDists = nan(8,length(opt.subjs_all),5);
medDists = nan(8,5);
SEM = nan(8,5);
for aa = 1:5
    for bb = 1:length(opt.subjs_all)
        indivDists(1,bb,aa) = median(dists{aa}{bb}.rp);
        indivDists(2,bb,aa) = median(dists{aa}{bb}.rple_eeg);
        for cc = 1:6
            indivDists(cc+2,bb,aa) = median(dists{aa}{bb}.rple_noise{cc});
        end
    end
    medDists(:,aa) = nanmedian(indivDists(:,:,aa),2);
    for cc = 1:8
        SEM(cc,aa) = 1.2533*(std(indivDists(cc,~isnan(indivDists(cc,:,aa)),aa))/...
            sqrt(length(opt.subjs_all)-sum(isnan(indivDists(cc,:,aa)))));
    end
end

figure
hold on
h = cell(1,2);
h{2} = bar([1,2,3,4,5,6],[medDists(1,1);nan;nan;nan;nan;nan],...
    'LineWidth',1.5);
h{2}(1).BarWidth = .3;
h{2}(1).FaceColor = [1,1,1];
errorbar([medDists(1,1);nan;nan;nan;nan;nan],...
    [SEM(1,1),nan,nan,nan,nan,nan],'k','linestyle','none',...
    'HandleVisibility','Off','LineWidth',3,'CapSize',10);
h{1} = bar([1,2,3,4,5,6],[[nan;nan;nan;nan;nan;nan;nan],...
    medDists(2:8,1),medDists(2:8,2),medDists(2:8,3),medDists(2:8,4),...
    medDists(2:8,5)],'LineWidth',1.5);
for aa = 1:7
    h{1}(aa).FaceColor = Colours(aa,:);
    h{1}(aa).Parent.XColor = [0,0,0];
    h{1}(aa).Parent.YColor = [0,0,0];
end

for bb = 1:7
    barloc = getErrorBarCentre(6,7,bb);
    errorbar(barloc,[nan;medDists(bb+1,1);medDists(bb+1,2);...
        medDists(bb+1,3);medDists(bb+1,4);medDists(bb+1,5)],...
        [nan;SEM(bb+1,1);SEM(bb+1,2);SEM(bb+1,3);SEM(bb+1,4);...
        SEM(bb+1,5)],'k','linestyle','none',...
        'HandleVisibility','Off','LineWidth',3,'CapSize',10);
end
xticks([1,2,3,4,5,6]);
xticklabels({'','1','2','3','4','5'});
lab = xlabel('Threshold');
lab.Position(1) = 4;
ylabel({'Similarity Score';'(Euclidean Distance to the RP Template)'});
leg = {'RPs'};
legEntries = cat(2,leg,legEntries);
legend(legEntries,'Location','northwest');
h{1}(1).Parent.FontSize = 20;
h{1}(1).Parent.FontWeight = 'bold';

% checks for significant differences between the similarity scores of RPs
% and RPLEs
sigdiff = cell(1,2);
sigdiff{1} = nan(5,6);
sigdiff{2} = nan(5,6);
for aa = 1:5
    for cc = 1:6
        [~,sigdiff{1}(aa,cc)] = kstest2(indivDists(2,:,aa),...
            indivDists(cc+2,:,aa));
        sigdiff{2}(aa,cc) = ranksum(indivDists(2,:,aa),...
            indivDists(cc+2,:,aa));
    end
end
tables.tables{2} = array2table(sigdiff{1},'VariableNames',varNames);
tables.tables{3} = array2table(sigdiff{2},'VariableNames',varNames);
tables.desc{2} = 'kstest2: median Euclidean distances (to the RP template) of the EEG RPLEs and the noise RPLEs';
tables.desc{3} = 'ranksum: median Euclidean distances (to the RP template) of the EEG RPLEs and the noise RPLEs';


%% Temporal plots

ColourOrd = [0,0,0;.4,.4,.4;.47,.67,.19;.64,.08,.18];
h = cell(1,5);
Ys = nan(5,2);
figure
for aa = 1:5
    for bb = 1:length(opt.subjs_all)
        sets{aa}.combepos{1}{bb} = ...
            proc_average(sets{aa}.combepos{1}{bb},'Stats',1);
    end
    epos = cell(1,8);
    for cc = 1:8
        epos{cc} = cell(1,length(opt.subjs_all));
        toRem = [];
        for bb = 1:length(opt.subjs_all)
            epos{cc}{bb} = proc_selectClasses(sets{aa}.combepos{1}{bb},...
                sets{aa}.combepos{1}{bb}.className{cc});
            if isnan(sum(epos{cc}{bb}.x))
                toRem = cat(2,toRem,bb);
            end
        end
        epos{cc}(toRem) = [];
        epos{cc} = proc_grandAverage(epos{cc},'Average','NWeighted');
    end
    gaepo = epos{1};
    gaepo.y = eye(8);
    gaepo.className = sets{1}.combepos{1}{1}.className;
    gaepo.className{1} = 'movement onset';
    gaepo.className{2} = 'rple eeg';
    RPstats = gaepo.se(:,:,2);
    gaepo.se(:,:,2) = gaepo.se(:,:,1);
    gaepo.se(:,:,1) = RPstats;
    RPstats = gaepo.p(:,:,2);
    gaepo.p(:,:,2) = gaepo.p(:,:,1);
    gaepo.p(:,:,1) = RPstats;
    RPstats = gaepo.sgnlogp(:,:,2);
    gaepo.sgnlogp(:,:,2) = gaepo.sgnlogp(:,:,1);
    gaepo.sgnlogp(:,:,1) = RPstats;
    for cc = 2:8
        gaepo.x = cat(3,gaepo.x,epos{cc}.x);
    end
    gaepo = proc_selectClasses(gaepo,{'rple eeg','movement onset',...
        'rple noise Per','rple noise 175'});
    
    subplot(3,2,aa)
    h{aa} = plot_channel(gaepo,'Cz','ColorOrder',ColourOrd,...
        'RefCol',[1,1,1]);
    h{aa}.ax.FontSize = 15;
    h{aa}.ax.FontWeight = 'bold';
    h{aa}.leg.Location = 'southwest';
    h{aa}.leg.FontWeight = 'bold';
    h{aa}.leg.FontSize = 12;
    h{aa}.ax.XColor = [0,0,0];
    h{aa}.ax.YColor = [0,0,0];
    h{aa}.XLabel.String = 'Time (ms)';
    h{aa}.XLabel.FontSize = 15;
    h{aa}.YLabel.String = 'Amplitude (\muV)';
    h{aa}.YLabel.FontSize = 15;
    h{aa}.XLabel.FontWeight = 'bold';
    h{aa}.YLabel.FontWeight = 'bold';
    for bb = 1:4
        h{aa}.plot(bb).LineWidth = 3;
    end
    h{aa}.title.String = '';
    Ys(aa,:) = h{aa}.ax.YLim;
end
Ys(1,1) = min(Ys(:,1));
Ys(1,2) = max(Ys(:,2));
Ys = Ys(1,:);
for aa = 1:5
    h{aa}.ax.YLim = Ys;
    h{aa}.leg.String = {'RP','EEG','Personal \alpha','\alpha = 1.75'};
end


end
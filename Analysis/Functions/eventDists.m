function [dists,tables] = eventDists(Cout_p1,Cout_rt,sets,...
    thresholds)
% eventDists - Gets the Euclidean distances between the RP template and 
%     the: single-trial RPs; single-event RPLEs (phase 1 and phase RT). 
%     Does so for each mode and each threshold. Plots the similarity scores
% 
% Argument(s):
%  Cout_p1 -  phase 1 classifier output for all subjects and all modes
%  Cout_rt -  phase RT classifier output for all subjects and all modes
%  sets -     contains the epoched RPs and RPLEs as well as information 
%             such as the number of RPLEs in the data for all subjects, for 
%             each mode and each threshold
%  thresholds - thresholds for each mode and subject used to identify RPLEs
% 
% Returns:
%  dists -    similarity scores of the RPs and RPLEs for each subject, each 
%             mode, and each threshold
%  tables -   tables of results for the statistical tests (kstest2, 
%             signtest) which tested whether the similarity scores of the 
%             RPLEs were significantly different to those of the RPs
% 
% Author(s): Thomas Binns, 2020


global opt

dists = cell(1,length(sets));
for aa = 1:length(sets)
    selthresh = cell(1,3);
    for bb = 1:3
        selthresh{bb} = thresholds{bb}(aa,:);
    end
    rpleepos = modeAnalysis(Cout_p1,Cout_rt,'thresh',selthresh,...
        'noAvg',0);
    % gets RP template for each subject
    template = cell(1,3);
    avgepos = cell(1,3);
    for bb = 1:3
        template{bb} = cell(1,length(opt.subjs_all));
        avgepos{bb} = cell(1,length(opt.subjs_all));
        for cc = 1:length(opt.subjs_all)
            avgepos{bb}{cc} = proc_average(rpleepos{bb}{cc});
            template{bb}{cc} = avgepos{bb}{cc}.x(:,:,...
                find(strcmp(avgepos{bb}{cc}.className,'movement onset')));
        end
    end
    dists{aa} = eucDist(rpleepos,template);
end

% gets medians of Euclidean distances for each participant
meddists.rp = nan(length(opt.subjs_all),3,5);
meddists.rple_p1 = meddists.rp;
meddists.rple_rt = meddists.rp;
for aa = 1:5
    for bb = 1:3
        for cc = 1:length(opt.subjs_all)
            meddists.rp(cc,bb,aa) = median(dists{aa}{bb}{cc}.rp);
            meddists.rple_p1(cc,bb,aa) = median(dists{aa}{bb}{cc}.rple_p1);
            meddists.rple_rt(cc,bb,aa) = median(dists{aa}{bb}{cc}.rple_rt);
        end
    end
end

% gets and plots grand medians
allmeds = nan(3,3,5);
for aa = 1:5
    for bb = 1:3
        allmeds(1,bb,aa) = median(meddists.rp(:,bb,aa));
        allmeds(2,bb,aa) = median(meddists.rple_p1(:,bb,aa));
        allmeds(3,bb,aa) = median(meddists.rple_rt(:,bb,aa));
    end
end

% gets standard error of medians
SEM = nan(size(allmeds));
for aa = 1:5
    for bb = 1:3
        SEM(1,bb,aa) = 1.2533*(std(meddists.rp(:,bb,aa))/...
            sqrt(length(meddists.rp(:,bb,aa))));
        SEM(2,bb,aa) = 1.2533*(std(meddists.rple_p1(:,bb,aa))/...
            sqrt(length(meddists.rple_p1(:,bb,aa))));
        SEM(3,bb,aa) = 1.2533*(std(meddists.rple_rt(:,bb,aa))/...
            sqrt(length(meddists.rple_rt(:,bb,aa))));
    end
end

h = cell(2,3);
ttl = {'Spatio-temporal Mode','Temporal Mode','Spatial Mode'};
figure
for aa = 1:3
    subplot(1,3,aa)
    hold on
    h{2,aa} = bar([1,2,3,4,5,6],[allmeds(1,aa,1);nan;nan;nan;nan;nan],...
        'LineWidth',1.5);
    h{2,aa}(1).FaceColor = [1,1,1];
    errorbar([allmeds(1,aa,1);nan;nan;nan;nan;nan],...
        [SEM(1,aa,1),nan,nan,nan,nan,nan],'k','linestyle','none',...
        'HandleVisibility','Off','LineWidth',3,'CapSize',10);
    h{1,aa} = bar([1,2,3,4,5,6],[[nan,nan];allmeds(2:3,aa,1)';...
        allmeds(2:3,aa,2)';allmeds(2:3,aa,3)';allmeds(2:3,aa,4)';...
        allmeds(2:3,aa,5)'],'LineWidth',1.5);
    h{1,aa}(1).FaceColor = [.5,.5,.5];
    h{1,aa}(2).FaceColor = [.8,.8,.8];
    groupwidth = 2/3.5;
    for bb = 1:2
        barloc = (1:6)-groupwidth/2+(2*bb-1)*groupwidth/(2*2);
        errorbar(barloc,[nan;allmeds(bb+1,aa,1);allmeds(bb+1,aa,2);...
            allmeds(bb+1,aa,3);allmeds(bb+1,aa,4);allmeds(bb+1,aa,5)],...
            [nan,SEM(bb+1,aa,1),SEM(bb+1,aa,2),SEM(bb+1,aa,3),...
            SEM(bb+1,aa,4),SEM(bb+1,aa,5)],'k','linestyle','none',...
            'HandleVisibility','Off','LineWidth',3,'CapSize',10);
    end
    xticks([1,2,3,4,5,6]);
    xticklabels({'','1','2','3','4','5'});
    title(ttl{aa});
    lab = xlabel('Threshold');
    lab.Position(1) = 4;
    if aa == 1
        ylabel('Euclidean Distance to the RP Template');
    end
    if aa == 2
        legend('RPs','Phase 1 RPLEs','Phase RT RPLEs','Location',...
            'northwest');
    end
    h{1,aa}(1).Parent.FontSize = 20;
    h{1,aa}(1).Parent.FontWeight = 'bold';
    h{1,aa}(1).Parent.XColor = [0,0,0];
    h{1,aa}(1).Parent.YColor = [0,0,0];
end

% checks for significant differences between the Euclidean distances of
% RPLEs and RPs
sigdiff = nan(5,6,2);
for aa = 1:5
    for bb = 1:3
        [~,sigdiff(aa,bb,1)] = kstest2(meddists.rp(:,bb,aa),...
            meddists.rple_p1(:,bb,aa));
        [~,sigdiff(aa,bb+3,1)] = kstest2(meddists.rp(:,bb,aa),...
            meddists.rple_rt(:,bb,aa));
        sigdiff(aa,bb,2) = signtest(meddists.rp(:,bb,aa),...
            meddists.rple_p1(:,bb,aa));
        sigdiff(aa,bb+3,2) = signtest(meddists.rp(:,bb,aa),...
            meddists.rple_rt(:,bb,aa));

    end
end
tables.data = cell(1,2);
tables.info = cell(1,2);
tables.info{1} = 'kstest2: median rp distances & median rple_p1 distances, median rple_rt distances';
tables.info{2} = 'signtest: median rp distances & median rple_p1 distances, median rple_rt distances';
vars = {'P1 ST','P1 T','P1 S','RT ST','RT T','RT S'};
rows = {'1','2','3','4','5'};
tables.data{1} = array2table(sigdiff(:,:,1),'VariableNames',vars,...
    'RowNames',rows);
tables.data{2} = array2table(sigdiff(:,:,2),'VariableNames',vars,...
    'RowNames',rows);


end
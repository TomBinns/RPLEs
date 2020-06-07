function plotEEGAnalysis(sets)
% plotEEGAnalysis - Plots results for the EEG data (number of RPLEs,
%     appearance of RPLEs)
%     
% Argument(s):
%  sets -     contains the epoched RPs and RPLEs as well as information 
%             such as the number of RPLEs in the data for all subjects, 
%             for each mode and each threshold
%  
% Author(s): Thomas Binns, 2020 

global opt

%% plots number of events per second

h = cell(1,3);
ttl = {'Spatio-temporal Mode','Temporal Mode','Spatial Mode'};
figure
for aa = 1:3
    SEM = cell(1,5);
    for bb = 1:5
        SEM{bb}(1,:) = std(sets{bb}.Eps{aa}(1,:))/...
            sqrt(size(sets{bb}.Eps{aa},2));
        SEM{bb}(2,:) = std(sets{bb}.Eps{aa}(2,:))/...
            sqrt(size(sets{bb}.Eps{aa},2));
    end
    subplot(1,3,aa)
    hold on
    h{aa} = bar([1,2,3,4,5],[mean(sets{1}.Eps{aa},2),...
        mean(sets{2}.Eps{aa},2),mean(sets{3}.Eps{aa},2),...
        mean(sets{4}.Eps{aa},2),mean(sets{5}.Eps{aa},2)],'LineWidth',...
        1.5);
    h{aa}(1).Parent.FontSize = 20;
    h{aa}(1).Parent.FontWeight = 'bold';
    for bb = 1:2
        centre = (1:5)-(2/3.5)/2+(2*bb-1)*(2/3.5)/4;
        errorbar(centre,[mean(sets{1}.Eps{aa}(bb,:)),...
            mean(sets{2}.Eps{aa}(bb,:)),mean(sets{3}.Eps{aa}(bb,:)),...
            mean(sets{4}.Eps{aa}(bb,:)),mean(sets{5}.Eps{aa}(bb,:))],...
            [SEM{1}(bb),SEM{2}(bb),SEM{3}(bb),SEM{4}(bb),SEM{5}(bb)],...
            'k','linestyle','none','LineWidth',3,'CapSize',10);
    end
    xticks([1,2,3,4,5]);
    xticklabels({'1','2','3','4','5'});
    title(ttl{aa});
    xlabel('Threshold')
    if aa == 1
        ylabel('RPLEs/s');
        legend('Phase 1 RPLEs','Phase RT RPLEs','Location','northeast');
    end
end
sgtitle('Number of RPLEs per second for different thresholds, phases and modes');


%% Spatio-temporal plots
[~,~,mnt] = loadData(opt.subjs_all{1},'Phase1');

h = cell(1,5);
Ys = nan(5,2);
for aa = 1:5
    gaepos{1} = proc_grandAverage(sets{aa}.rpleepos{1},'Average',...
        'Nweighted');
    figure
    h{aa} = grid_plot(gaepos{1},mnt,'ShrinkAxes',[0.9,0.9],'PlotStat',...
        'sem');
    Ys(aa,:) = h{aa}.scale.ax.YLim;
    close;
end
Ys(1,1) = min(Ys(:,1));
Ys(1,2) = max(Ys(:,2));
Ys = Ys(1,:);
for aa = 1:5
    gaepos{1} = proc_grandAverage(sets{aa}.rpleepos{1},'Average',...
        'Nweighted');
    gaepos{1}.y([1,2],:) = gaepos{1}.y([2,1],:);
    gaepos{1}.className = {'RP','Phase 1 RPLE','Phase RT RPLE'};
    figure
    h{aa} = grid_plot(gaepos{1},mnt,'ShrinkAxes',[0.83,0.83],...
        'ScalePolicy',Ys,'ColorOrder',[0,0,0;.4,.4,.4;.7,.7,.7],...
        'RefCol',[1,1,1],'AxisTitleFontWeight','bold',...
        'AxisTitleFontSize',15);
    h{aa}.leg.FontSize = 15;
    h{aa}.leg.FontWeight = 'bold';
end


%% Temporal plots
h = cell(1,5);
Ys = nan(5,2);
figure
for aa = 1:5
    gaepos{2} = proc_grandAverage(sets{aa}.rpleepos{2},'Average',...
        'Nweighted');
    gaepos{2}.y([1,2],:) = gaepos{2}.y([2,1],:);
    gaepos{2}.className = {'RP','Phase 1 RPLE','Phase RT RPLE'};
    subplot(3,2,aa)
    h{aa} = plot_channel(gaepos{2},'1',...
        'ColorOrder',[0,0,0;.4,.4,.4;.7,.7,.7],'RefCol',[1,1,1]);
    h{aa}.ax.FontSize = 15;
    h{aa}.ax.FontWeight = 'bold';
    h{aa}.leg.Location = 'southwest';
    h{aa}.leg.FontWeight = 'bold';
    h{aa}.leg.FontSize = 12;
    h{aa}.XLabel.String = 'Time (ms)';
    h{aa}.XLabel.FontWeight = 'bold';
    h{aa}.XLabel.FontSize = 15;
    h{aa}.YLabel.String = 'Amplitude (\muV)';
    h{aa}.YLabel.FontWeight = 'bold';
    h{aa}.YLabel.FontSize = 15;
    for bb = 1:3
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
end


%% Spatial plots

h = cell(1,5);
for aa = 1:5
    h{aa} = cell(1,3);
end
Ys = nan(5,2);

for aa = 1:5
    gaepos{3} = proc_grandAverage(sets{aa}.rpleepos{3},'Average',...
        'Nweighted');
%     garsq = proc_grandAverage(sets{aa}.rsqs,'Stats',1);
%     garsq.x(garsq.p>0.05) = 0;
    Ys(aa,1) = max(abs(gaepos{3}.x(:)))*-1;
%     Ys(aa,2) = max(abs(garsq.x(:)))*-1;
%     Ys(aa,1) = min(gaepos{3}.x(:));
%     Ys(aa,2) = max(gaepos{3}.x(:));
end
Ys(1,1) = min(Ys(:,1));
Ys(2,1) = min(Ys(:,2));
Ys(1,2) = Ys(1,1)*-1;
Ys(2,2) = Ys(2,1)*-1;
Ys = Ys(1:2,:);

mnt = mnt_adaptMontage(mnt,gaepos{1});
for aa = 1:5
    gaepos{3} = proc_grandAverage(sets{aa}.rpleepos{3},'Average',...
        'Nweighted');
%     gaepos{3}.y([1,2],:) = gaepos{3}.y([2,1],:);
    gaepos{3}.className = {'Phase 1 RPLE','Phase RT RPLE','RP'};
    figure
    for bb = 1:3
        subplot(1,3,bb)
        hold on
        h{aa}{bb} = plot_scalp(mnt,gaepos{3}.x(end,:,bb),defopt_scalp_r...
            ('Extrapolation',1,'ExtrapolateToZero',1),'CLim',Ys(1,:));
        title(gaepos{3}.className{bb},'FontSize',17);
        if bb < 3
            h{aa}{bb}.cb.Visible = 'off';
        else
            h{aa}{bb}.cb.Color = [0,0,0];
            h{aa}{bb}.cb.FontSize = 15;
            h{aa}{bb}.cb.FontWeight = 'bold';
            h{aa}{bb}.cb.TicksMode = 'auto';
            h{aa}{bb}.cb.Position = [0.93,0.31,0.01,0.41];
            h{aa}{bb}.cb.Label.String = 'Amplitude (\muV)';
            h{aa}{bb}.cb.Label.FontSize = 17;
            h{aa}{bb}.cb.Label.Color = [0,0,0];
        end
    end
end


end
function [EpsTable,distsTables] = analysis_EEG
% analysis_EEG - Find RPLEs in the EEG data, gets information on the number
%     of events and the similarity scores of the events
% 
% Returns:
%  EpsTable - table of results for the ttest for each mode and threshold 
%             which tested whether there were more events in phase 1 than 
%             phase RT
%  distsTables - tables of results for the kstest2 and signtest for each 
%             mode and threshold which tested whether the similarity scores
%             of the RPLEs (phase 1; phase RT) were different to the
%             similarity scores of the RPs
% 
% Author(s): Thomas Binns, 2020

%% sets up analysis

global opt
opt.subjs_all = {'VPtae','VPtaf','VPtah','VPtai',...
                 'VPtal','VPtam','VPtan','VPtao','VPtap','VPtaq',...
                 'VPtar','VPtas','VPtat','VPtau','VPtay',...
                 'VPtaz','VPtba'};

opt.clab_load = {'F1','Fz','F2',...
             'FC3','FC2','FC1','FC4',...
             'C3','C1','Cz','C2','C4',...
             'CP3','CP1','CPz','CP2','CP4',...
             'P1','Pz','P2'};
opt.baseln_len = 100;
opt.baseln_pos = 'beginning';
opt.untilTime = -1000;
opt.WindowStartPoint = -300;
opt.WindowEndPoint = 200;
opt.method = 'setval';
opt.beta = 1;
opt.Cival = [1000,NaN];
opt.epochSegment = [-1000,0];


%% loads Cout

CoutDir = append(opt.dataDir,'\Cout');
loadname = fullfile(CoutDir,'Cout_p1.mat');
load(loadname);
loadname = fullfile(CoutDir,'Cout_rt.mat');
load(loadname);


%% runs 5-threshold analysis

% loads threshold
threshDir = append(opt.dataDir,'\Threshold\thresholds.mat');
load(threshDir);
%thresholds = findThresholds(Cout_p1);

sets = cell(1,5);

for aa = 1:5
    selthresh = cell(1,3);
    for bb = 1:3
        selthresh{bb} = thresholds{bb}(aa,:);
    end
    [rpleepos,Eps,rpleTs,rsqs] = modeAnalysis(Cout_p1,...
        Cout_rt,'thresh',selthresh);
    sets{aa}.rpleepos = rpleepos;
    sets{aa}.Eps = Eps;
    sets{aa}.rpleTs = rpleTs;
    sets{aa}.rsqs = rsqs;
end

EpsTab = table(nan(5,1),nan(5,1),nan(5,1),'VariableNames',...
    {'ST','T','S'},'RowNames',{'Threshold 1','Threshold 2',...
    'Threshold 3','Threshold 4','Threshold 5'});
for aa = 1:5
    for bb = 1:3
        [~,EpsTab{aa,bb}] = ttest(sets{aa}.Eps{bb}(1,:),...
            sets{aa}.Eps{bb}(2,:),'Tail','right');
    end
end

EpsTable.data = EpsTab;
EpsTable.info = 'ttest: NRPLEs phase 1 > NRPLEs phase RT';

plotEEGAnalysis(sets);

[dists,distsTables] = eventDists(Cout_p1,Cout_rt,sets,thresholds,0);  

end

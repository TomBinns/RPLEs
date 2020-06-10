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
                 'VPtaz','VPtba'}; % subject list
opt.clab_load = {'F1','Fz','F2',...
             'FC3','FC2','FC1','FC4',...
             'C3','C1','Cz','C2','C4',...
             'CP3','CP1','CPz','CP2','CP4',...
             'P1','Pz','P2'}; % channels to load
opt.baseln_len = 100; % length of baseline (ms)
opt.baseln_pos = 'beginning'; % baseline gathered from first baseln_len ms
opt.untilTime = -1000; % start of the RP period / time until which 
    % threshold crossings count as RPLEs in each trial (in ms; e.g. -1000 =
    % 1000 ms before movement onset)
opt.Cival = [1000,NaN]; % Cival(1): interval that must occur between
    % threshold crossings (in ms); threshold crossings that occur within 
    % this time of a prior crossing are removed
    % Cival(2): start of the RP 'hit' window; thresholds crossings that
    % occur within (1) prior to this period are removed; if NaN, crossings
    % not removed within (1) of the 'hit' window
opt.epochSegment = [-1000,0]; % segment of data to epoch, from (1) ms prior
    % to the event marker to (2) ms after the event marker

% not exactly necessary, just as a precautiong
opt.WindowStartPoint = -300; % start of the RP 'hit' window (in ms)
opt.WindowEndPoint = 200; % end of the RP 'hit' window (in ms)
opt.beta = 1; % beta value used in the F-measure (see Powers, 2011)
% 'hit' window refers to the time period during which a threshold crossing
% counts as a RP


%% loads classifier outputs

CoutDir = append(opt.dataDir,'\Cout');
loadname = fullfile(CoutDir,'Cout_p1.mat');
load(loadname);
loadname = fullfile(CoutDir,'Cout_rt.mat');
load(loadname);


%% runs analysis to find RPLEs in the EEG data for all 3 modes and all 5 
% thresholds

% generates thresholds
% thresholds = findThresholds(Cout_p1);

% loads threshold
threshDir = append(opt.dataDir,'\Threshold\thresholds.mat');
load(threshDir);

% runs analysis on EEG data
sets = cell(1,5);
for aa = 1:5
    selthresh = cell(1,3);
    for bb = 1:3
        selthresh{bb} = thresholds{bb}(aa,:);
    end
    [rpleepos,Eps] = modeAnalysis(Cout_p1,...
        Cout_rt,'thresh',selthresh);
    sets{aa}.rpleepos = rpleepos;
    sets{aa}.Eps = Eps;
end

% checks if there are significantly more events in phase 1 than phase RT
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

% plots results for number of events and appearence of events
plotEEGAnalysis(sets);

% plots results for similarity scores of events and checks if scores of
% RPLEs are significantly different to scores of RPs
[~,distsTables] = eventDists(Cout_p1,Cout_rt,sets,thresholds);  

end

function [tables] = analysis_noise
% analysis_noise - Find RPLEs in the noise data, gets information on the 
%     number of events and the similarity scores of the events
% 
% Returns:
%  tables -   results tables for the statistical tests:
%              ttest2 - compares number of RPLEs in the EEG data to the
%                       noise data
%              kstest2 and ranksum - compare the similarity scores of EEG 
%                                    RPLEs to the noise RPLEs
% 
% Author(s): Thomas Binns, 2020


global opt

opt.subjs_all = {'VPtae','VPtaf','VPtah','VPtai',...
                 'VPtal','VPtam','VPtan','VPtao','VPtap','VPtaq',...
                 'VPtar','VPtas','VPtat','VPtau','VPtay',...
                 'VPtaz','VPtba'}; % subject list
opt.clab_load = {'Cz'}; % channels to load
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


%% loads classifier output

CoutDir = fullfile(opt.dataDir,'Cout');
loadname = fullfile(CoutDir,'Cout_Cz_p1.mat');
load(loadname);
loadname = fullfile(CoutDir,'Cout_Cz_rt.mat');
load(loadname);


%% runs analysis to find RPLEs in EEG and noise data using the Cz-only mode
% of analysis and all 5 thresholds

% generates thresholds
% thresholds = findCzThresholds(Cout_p1);

% loads thresholds
threshDir = fullfile(opt.dataDir,'Threshold\CzThresholds.mat');
load(threshDir);

% runs analysis on EEG data
sets = cell(1,5);
for aa = 1:5
    selthresh = cell(1,3);
    for bb = 1:1
        selthresh{1} = thresholds(aa,:);
    end
    [rpleepos,Eps,~] = modeAnalysisCz(Cout_p1,Cout_rt,...
        'thresh',selthresh);
    sets{aa}.rpleepos = rpleepos;
    sets{aa}.Eps = Eps;
    sets{aa}.rpleTs = rpleTs;
end

% runs analysis on noise data
CoutTypes = {'Cz_000','Cz_050','Cz_100','Cz_150','Cz_Personal','Cz_175'};
cntTypes = {'Noise000','Noise050','Noise100','Noise150','NoisePer',...
    'Noise175'};
for bb = 1:6
    % loads noise data
    loadname = sprintf('Cout_%s.mat',CoutTypes{bb});
    Coutname = fullfile(CoutDir,loadname);
    load(Coutname);
    for aa = 1:5
        [sets{aa}.noiseepos(bb),sets{aa}.noiseEps(bb)] = ...
            modeAnalysisNoise(Cout_noise,thresholds(aa,:),...
            cntTypes{bb});
    end
end

% plots analysis (number of events, similarity scores, appearence of 
% events) and performs statistical analysis (number of events, similarity
% scores)
[~,tables] = plotNoiseAnalysis(sets,thresholds);


end
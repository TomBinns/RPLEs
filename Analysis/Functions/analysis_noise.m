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
                 'VPtaz','VPtba'};

opt.clab_load = {'Cz'};
opt.baseln_len = 100;
opt.baseln_pos = 'beginning';
opt.untilTime = -1000;
opt.WindowStartPoint = -300;
opt.WindowEndPoint = 200;
opt.beta = 1;
opt.Cival = [1000,NaN];
opt.epochSegment = [-1000,0];


%% loads Cout

CoutDir = fullfile(opt.dataDir,'Cout');
loadname = fullfile(CoutDir,'Cout_Cz_p1.mat');
load(loadname);
loadname = fullfile(CoutDir,'Cout_Cz_rt.mat');
load(loadname);


%% runs analysis

% thresholds = findCzThresholds(Cout_p1);
threshDir = fullfile(opt.dataDir,'Threshold\CzThresholds.mat');
load(threshDir);

sets = cell(1,5);
for aa = 1:5
    selthresh = cell(1,3);
    for bb = 1:1
        selthresh{1} = thresholds(aa,:);
    end
    [rpleepos,Eps,rpleTs] = modeAnalysisCz(Cout_p1,Cout_rt,...
        'thresh',selthresh);
    sets{aa}.rpleepos = rpleepos;
    sets{aa}.Eps = Eps;
    sets{aa}.rpleTs = rpleTs;
end

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

[~,tables] = plotNoiseAnalysis(sets,thresholds);

end
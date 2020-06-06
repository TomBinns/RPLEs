%% Gets classifier outputs for the simulated noise data

global opt

% specify directory where data is located
opt.dataDir = 'C:\Users\tomth\Documents\Coding\RPLE_data';

opt.subjs_all = {'VPtae','VPtaf','VPtah','VPtai',...
    'VPtal','VPtam','VPtan','VPtao','VPtap','VPtaq',...
    'VPtar','VPtas','VPtat','VPtau','VPtay',...
    'VPtaz','VPtba'};
Ns = length(opt.subjs_all);

opt.clab_load = {'F1','Fz','F2',...
             'FC3','FC2','FC1','FC4',...
             'C3','C1','Cz','C2','C4',...
             'CP3','CP1','CPz','CP2','CP4',...
             'P1','Pz','P2'};


%%

ival_fv = [-1000 0];
bsln_len = 100;
bsln_pos = 'beginning';
ivals_jm = [-1000  -900;
            -900   -800;
            -800   -700;
            -700   -600;
            -600   -500;
            -500   -400;
            -400   -300;
            -300   -200;
            -200   -100;
            -100      0];

% temporal features
proc{1}.train = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_meanAcrossChannels}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{1}.apply = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{1}.ival = ival_fv;


phaseName = ''; % noise data you want to get classifier output of
if ~strcmp(phaseName,'000') && ~strcmp(phaseName,'050') &&...
        ~strcmp(phaseName,'100') && ~strcmp(phaseName,'150') &&...
        ~strcmp(phaseName,'Per') && ~strcmp(phaseName,'175')
    error('unrecognised alpha value')
end

Np = length(proc);
Cout_noise = cell(Ns,Np);
for ii = 1:Ns
    
    % load phase 1 data
    [cnt1,mrk1] = loadData(opt.subjs_all{ii},'Phase1');
    cnt1 = proc_selectChannels(cnt1,{'Cz'});
    trl = getTrialMarkers(mrk1,'movement onset');
    mrk1 = mrk_selectEvents(mrk1,[trl{:}]);
    mrk1 = mrk_selectClasses(mrk1,{'trial start','movement onset','trial end'});
    
    % load noise data
    [cnt2,mrk2] = loadData(opt.subjs_all{ii},phaseName);
    
    for jj = 1:Np
        
        mrk_ = mrk_selectClasses(mrk1,{'trial start','movement onset'});
        fv = proc_segmentation(cnt1,mrk_,proc{jj}.ival);
        
        % train classifier for sliding noise cout
        fv = proc_chain(fv,proc{jj}.train);
        C = train_RLDAshrink(fv.x,fv.y);
        
        % sliding cout
        Cout_noise{ii,jj} = proc_slidingClassi(cnt2,mrk2,proc{jj},C);
        
    end
    
end


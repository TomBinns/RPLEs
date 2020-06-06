%% Gets the accuracies of the classifiers for all modes, based on the loss
% (as determined by the ROC curves; accuracy (%) = (1-loss)*100)

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
ival_amp = [-400 0];

%%

% spatio-temporal features
proc{1}.train = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{1}.apply = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{1}.ival = ival_fv;

% temporal features
proc{2}.train = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_meanAcrossChannels}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{2}.apply = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_meanAcrossChannels}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{2}.ival = ival_fv;

% spatial features
proc{3}.train = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_meanAcrossTime,ival_amp}
                 {@proc_flaten}};
proc{3}.apply = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_meanAcrossTime,ival_amp}
                 {@proc_flaten}};
proc{3}.ival = ival_fv;

% Cz-only temporal features
proc{4}.train = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_selectChannels,{'Cz'}}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{4}.apply = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_selectChannels,{'Cz'}}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{4}.ival = ival_fv;


Np = length(proc);
lossROC = nan(Ns,Np);
for ii = 1:Ns
    
    % load phase 1 data
    [mrk1,cnt1] = loadData(opt.subjs_all{ii},'Phase1');
    cnt1 = proc_selectChannels(cnt1,opt.clab_load);
    trl = getTrialMarkers(mrk1,'movement onset');
    mrk1 = mrk_selectEvents(mrk1,[trl{:}]);
    mrk1 = mrk_selectClasses(mrk1,{'trial start','movement onset',...
        'trial end'});
    
    for jj = 1:Np
        
        mrk_ = mrk_selectClasses(mrk1,{'trial start','movement onset'});
        fv = proc_segmentation(cnt1,mrk_,proc{jj}.ival);
        
        % get cross-validated accuracies
        [~,~,Cout] = crossvalidation(fv,@train_RLDAshrink,...
            'SampleFcn',{@sample_leaveOneOut},'Proc',proc{jj});
        loss = loss_rocArea(fv.y,Cout);
        lossROC(ii,jj) = (1-loss)*100;
        
    end
    
end
         
         
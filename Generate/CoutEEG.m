%% Gets classifier output for the EEG data

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

%% Gets classifier outputs for the EEG data for the spatio-temporal, 
% temporal, and spatial modes of analysis

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


rt = 0; % if rt = 0, trained classifier applied to phase 1 data 
        % if rt = 1, trained classifier applied to phase RT data

Np = length(proc);
Cout_p1 = cell(Ns,Np);
Cout_rt = cell(Ns,Np);
for ii = 1:Ns
    
    % load phase 1 data
    [mrk1,cnt1] = loadData(opt.subjs_all{ii},'Phase1');
    cnt1 = proc_selectChannels(cnt1,opt.clab_load);
    trl = getTrialMarkers(mrk1,'movement onset');
    mrk1 = mrk_selectEvents(mrk1,[trl{:}]);
    mrk1 = mrk_selectClasses(mrk1,{'trial start','movement onset',...
        'trial end'});
    
    if rt == 1
        % load rt data
        [mrk2,cnt2] = loadData(opt.subjs_all{ii},'RT');
        cnt2 = proc_selectChannels(cnt2,opt.clab_load);
        trl = getTrialMarkers(mrk2,'go signal');
        mrk2 = mrk_selectEvents(mrk2,[trl{:}]);
        mrk2 = mrk_selectClasses(mrk2,{'trial start','go signal',...
            'trial end'});
    end
    
    for jj = 1:Np
        
        mrk_ = mrk_selectClasses(mrk1,{'trial start','movement onset'});
        fv = proc_segmentation(cnt1,mrk_,proc{jj}.ival);
        
        if rt == 1
            % train classifier for sliding rt cout
            fv = proc_chain(fv,proc{jj}.train);
            C = train_RLDAshrink(fv.x,fv.y);
            % sliding cout
            Cout_rt{ii,jj} = proc_slidingClassi(cnt2,mrk2,proc{jj},C);
        else
            % sliding cout
            Cout_p1{ii,jj} = proc_slidingClassi(cnt1,mrk1,proc{jj});
        end
        
    end
    
end
         
         
%% Gets classifier outputs for the EEG data for the Cz-only mode of 
% analysis

% temporal features
proc{1}.train = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{1}.apply = {{@proc_baseline,bsln_len,bsln_pos}
                 {@proc_jumpingMeans,ivals_jm}
                 {@proc_flaten}};
proc{1}.ival = ival_fv;


rt = 1; % if rt = 0, trained classifier applied to phase 1 data 
        % if rt = 1, trained classifier applied to phase RT data

Np = length(proc);
Cout_p1 = cell(Ns,Np);
Cout_rt = cell(Ns,Np);
for ii = 1:Ns
    
    % load phase 1 data
    [mrk1,cnt1] = loadData(opt.subjs_all{ii},'Phase1');
    cnt1 = proc_selectChannels(cnt1,{'Cz'});
    trl = getTrialMarkers(mrk1,'movement onset');
    mrk1 = mrk_selectEvents(mrk1,[trl{:}]);
    mrk1 = mrk_selectClasses(mrk1,{'trial start','movement onset',...
        'trial end'});
    
    if rt == 1
        % load rt data
        [mrk2,cnt2] = loadData(opt.subjs_all{ii},'RT');
        cnt2 = proc_selectChannels(cnt2,{'Cz'});
        trl = getTrialMarkers(mrk2,'go signal');
        mrk2 = mrk_selectEvents(mrk2,[trl{:}]);
        mrk2 = mrk_selectClasses(mrk2,{'trial start','go signal',...
            'trial end'});
    end
    
    for jj = 1:Np
        
        mrk_ = mrk_selectClasses(mrk1,{'trial start','movement onset'});
        fv = proc_segmentation(cnt1,mrk_,proc{jj}.ival);
        
        if rt == 1
            % train classifier for rt cout
            fv = proc_chain(fv,proc{jj}.train);
            C = train_RLDAshrink(fv.x,fv.y);
            % sliding cout
            Cout_rt{ii,jj} = proc_slidingClassi(cnt2,mrk2,proc{jj},C);
        else
            % sliding cout
            Cout_p1{ii,jj} = proc_slidingClassi(cnt1,mrk1,proc{jj});
        end
        
    end
    
end


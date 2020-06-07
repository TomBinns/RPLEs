%% Gets classifier output for the EEG data

global opt

% specify directory where data is located
opt.dataDir = 'C:\Users\tomth\Documents\Coding\RPLE_data';

opt.subjs_all = {'VPtae','VPtaf','VPtah','VPtai',...
    'VPtal','VPtam','VPtan','VPtao','VPtap','VPtaq',...
    'VPtar','VPtas','VPtat','VPtau','VPtay',...
    'VPtaz','VPtba'}; % subject list
Ns = length(opt.subjs_all);
opt.clab_load = {'F1','Fz','F2',...
             'FC3','FC2','FC1','FC4',...
             'C3','C1','Cz','C2','C4',...
             'CP3','CP1','CPz','CP2','CP4',...
             'P1','Pz','P2'}; % channels to load
         
ival_fv = [-1000 0]; % periods to train classifier on
bsln_len = 100; % length of baseline (ms)
bsln_pos = 'beginning'; % what period to use for baseline
ivals_jm = [-1000  -900; % intervals for jumping means
            -900   -800;
            -800   -700;
            -700   -600;
            -600   -500;
            -500   -400;
            -400   -300;
            -300   -200;
            -200   -100;
            -100      0];
ival_amp = [-400 0]; % period for taking mean-across-time for spatial mode

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
    trl = getTrialMarkers(mrk1,'movement onset'); % only use data from
        % trials containing a movement
    mrk1 = mrk_selectEvents(mrk1,[trl{:}]);
    mrk1 = mrk_selectClasses(mrk1,{'trial start','movement onset',...
        'trial end'});
    
    if rt == 1
        % load phase RT data
        [mrk2,cnt2] = loadData(opt.subjs_all{ii},'RT');
        cnt2 = proc_selectChannels(cnt2,opt.clab_load);
        trl = getTrialMarkers(mrk2,'go signal'); % only use data from
            % trials containing a cue to move
        mrk2 = mrk_selectEvents(mrk2,[trl{:}]);
        mrk2 = mrk_selectClasses(mrk2,{'trial start','go signal',...
            'trial end'});
    end
    
    for jj = 1:Np
        
        mrk_ = mrk_selectClasses(mrk1,{'trial start','movement onset'});
        fv = proc_segmentation(cnt1,mrk_,proc{jj}.ival); % isolate data
            % preceding trial start (idle period) and movement onset
            % (RP period)
        
        if rt == 1
            % train classifier before getting phase RT classifier output
            fv = proc_chain(fv,proc{jj}.train);
            C = train_RLDAshrink(fv.x,fv.y);
            % get sliding classifier output
            Cout_rt{ii,jj} = proc_slidingClassi(cnt2,mrk2,proc{jj},C);
        else
            % get sliding classifier output
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
    cnt1 = proc_selectChannels(cnt1,{'Cz'}); % take only data from Cz
    trl = getTrialMarkers(mrk1,'movement onset'); % only use data from
        % trials containing a movement
    mrk1 = mrk_selectEvents(mrk1,[trl{:}]);
    mrk1 = mrk_selectClasses(mrk1,{'trial start','movement onset',...
        'trial end'});
    
    if rt == 1
        % load phase RT data
        [mrk2,cnt2] = loadData(opt.subjs_all{ii},'RT');
        cnt2 = proc_selectChannels(cnt2,{'Cz'}); % take only data from Cz
        trl = getTrialMarkers(mrk2,'go signal'); % only use data from
            % trials containing a cue to move
        mrk2 = mrk_selectEvents(mrk2,[trl{:}]);
        mrk2 = mrk_selectClasses(mrk2,{'trial start','go signal',...
            'trial end'});
    end
    
    for jj = 1:Np
        
        mrk_ = mrk_selectClasses(mrk1,{'trial start','movement onset'});
        fv = proc_segmentation(cnt1,mrk_,proc{jj}.ival); % isolate data
            % preceding trial start (idle period) and movement onset
            % (RP period)
        
        if rt == 1
            % train classifier for getting RT classifier output
            fv = proc_chain(fv,proc{jj}.train);
            C = train_RLDAshrink(fv.x,fv.y);
            % gets sliding classifier output
            Cout_rt{ii,jj} = proc_slidingClassi(cnt2,mrk2,proc{jj},C);
        else
            % gets sliding classifier output
            Cout_p1{ii,jj} = proc_slidingClassi(cnt1,mrk1,proc{jj});
        end
        
    end
    
end


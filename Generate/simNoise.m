%% Simulates single-channels' worth of noise with specific spectral
% densities (alpha values of the PSDs), scaled to match EEG data of
% channel Cz

% Author(s): Thomas Binns, 2020

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
opt.untilTime = -1000; % time until which RPLEs are found; only used if
    % generating trial markers for the noise


%% Gets alpha value of EEG data (mean of phase 1 and phase RT alpha values)

excludeFinal = 1000; % excludes RP period from phase 1 data
Hz = 30; % frequencies up to this value will be analysed
alphas = nan(Ns,2);
fits = nan(Ns,2,2);
yVals = nan(Ns,Hz,2);
yLim = [];
for aa = 1:Ns
    
    % phase 1
    data = [];
    
    % loads data
    [mrk,cnt] = loadData(opt.subjs_all{aa},'Phase1');
    cnt = proc_selectChannels(cnt,'Cz'); % uses only data from Cz
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset',...
        'trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,2000,{'trial start','movement onset'});
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset'});
    
    % isolates within-trial data (excluding the RP period)
    times = convTimes(mrk.time,100);
    idx = ones(size(times));
    idx(1:2:end) = 0;
    times(logical(idx)) = times(logical(idx)) - (excludeFinal/10);
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    cnt.x = data;
    
    % gets alpha value of subject's EEG data
    spec = proc_spectrum(cnt,[1,Hz],'Scaling','power');
    xvals = log10(spec.t);
    yvals = log10(spec.x);
    yvals = yvals';
    fit = polyfit(xvals,yvals,1); % log-log plots data and gets gradient
    alphas(aa,1) = fit(1)*-1;
    fits(aa,1,1) = fit(1);
    fits(aa,2,1) = fit(2);
    yVals(aa,:,1) = yvals;
    
    % phase RT
    data = [];
    
    % loads data
    [mrk,cnt] = loadData(opt.subjs_all{aa},'RT');
    cnt = proc_selectChannels(cnt,'Cz'); % uses only data from Cz
    mrk = mrk_selectClasses(mrk,{'trial start','go signal',...
        'trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,1000,{'trial start','go signal'});
    mrk = mrk_selectClasses(mrk,{'trial start','go signal'});
    
    % isolates within-trial data (excluding the RP period)
    times = convTimes(mrk.time,100);
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    cnt.x = data;
    
    % gets alpha value of subject's EEG data
    spec = proc_spectrum(cnt,[1,Hz],'Scaling','power');
    xvals = log10(spec.t);
    yvals = log10(spec.x);
    yvals = yvals';
    fit = polyfit(xvals,yvals,1); % log-log plots data and gets gradient
    alphas(aa,2) = fit(1)*-1;
    fits(aa,1,2) = fit(1);
    fits(aa,2,2) = fit(2);
    yVals(aa,:,2) = yvals;
    
end

% plots grand average alpha values of phase 1 and phase RT EEG data
meanalphas = mean(alphas,2);
meanYVals = mean(yVals,1);
meanFits = mean(fits,3);
meanValFits = polyval(mean(meanFits,1),xvals);
figure
hold on
plot((1:Hz),meanYVals(:,:,1),'LineWidth',5,'Color',[.5,.5,.5]);
plot((1:Hz),meanYVals(:,:,2),'LineWidth',5,'Color',[.8,.8,.8]);
h = plot((1:Hz),meanValFits,':','LineWidth',3,'Color',[0,0,0]);
h.Parent.FontSize = 20;
h.Parent.FontWeight = 'bold';
xlabel('Frequency (Hz)');
ylabel('Log_{10} Power (\muV^2/Hz)');
alphaLeg = sprintf('\\alpha \\approx %.2f',round(mean(meanalphas),2));
legend('Phase 1','Phase RT',alphaLeg);


%% Gets all WTs and randomly chooses times to create trial markers for
% noise data
% 
% % gets trial lengths of all trials and gets number of trials per subject
% % N.B. phase 1 data only
% Ntrials = nan(1,Ns);
% allWTs = cell(1,Ns);
% catWTs = [];
% for aa = 1:Ns
%     mrk = loadData(opt.subjs_all{aa},'Phase1');
%     mrk = mrk_selectClasses(mrk,{'trial start','movement onset',...
%         'trial end'});
%     trial_mrk = getTrialMarkers(mrk);
%     trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
%     mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
%     mrk = removeShortTrials(mrk,(1000+abs(opt.untilTime)),...
%         {'trial start','movement onset'});
%     allWTs{aa} = getWTs(mrk,{'trial start','movement onset'});
%     catWTs = cat(2,catWTs,allWTs{aa});
%     Ntrials(aa) = length(allWTs{aa});
% end
% catWTs = sort(catWTs);
% 
% % randomly chooses trial lengths
% idx = randi(length(catWTs),1,ceil(mean(Ntrials)));
% selWTs = catWTs(idx);
% 
% % creates new trial markers
% newmrk = mrk;
% newmrk.time = [];
% newmrk.y = repmat(eye(3),1,length(idx));
% zz = 1;
% for aa = 1:length(idx)
%     % trial start marker
%     if aa == 1
%         newmrk.time(zz) = 1;
%     else
%         newmrk.time(zz) = newmrk.time(zz-1) + 3000; % inter-trial period
%     end
%     % movement onset marker
%     newmrk.time(zz+1) = newmrk.time(zz) + selWTs(aa);
%     % trial end marker
%     newmrk.time(zz+2) = newmrk.time(zz+1) + 2000;
%     zz = zz+3;
% end
% mrk = newmrk;


%% Simulates noise data based on alpha values of EEG data

mrk = loadData(opt.subjs_all{1},'NoisePer');
newtimes = convTimes(mrk.time,100); % finds how long data needs to be

excludeFinal = 1000; % excludes RP period from phase 1 data

% gets standard deviation of channel Cz of EEG data (average std of phase 1
% and phase RT)
stds = nan(Ns,2);
for aa = 1:Ns
    
    % phase 1
    data = [];
    
    % loads data
    [mrk,cnt] = loadData(opt.subjs_all{aa},'Phase1');
    cnt = proc_selectChannels(cnt,'Cz'); % uses only data from Cz
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset',...
        'trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,2000,{'trial start','movement onset'});
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset'});
    
    % isolates within-trial data (excluding the RP period)
    times = convTimes(mrk.time,100);
    idx = ones(size(times));
    idx(1:2:end) = 0;
    times(logical(idx)) = times(logical(idx)) - (excludeFinal/10);
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    stds(aa,1) = std(data); % gets standard deviation of data
    
    % phase RT
    data = [];
    
    % loads data
    [mrk,cnt] = loadData(opt.subjs_all{aa},'RT');
    cnt = proc_selectChannels(cnt,'Cz'); % uses only data from Cz
    mrk = mrk_selectClasses(mrk,{'trial start','go signal',...
        'trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,1000,{'trial start','go signal'});
    mrk = mrk_selectClasses(mrk,{'trial start','go signal'});
    
    % isolates within-trial data (excluding the RP period)
    times = convTimes(mrk.time,100);
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    stds(aa,2) = std(data); % gets standard deviation of data
    
end
alphaVals = round(meanalphas,2);
stdVals = mean(stds,2);


saveDir = ''; % filepath to save simulated noise
for aa = 1:Ns
    % generates noise
    base = single(randn(newtimes(end),1));
    noiseOut = genNoise(alphaVals(aa),base);
    % scales noise
    scalednoise = noiseOut/std(noiseOut);
    scalednoise = scalednoise*stdVals(aa);
    % puts into cnt format and saves cnt
    cnt.x = scalednoise;
    cnt.clab = {'Cz'};
    name = sprintf('cnt_%s_Cz_Per',opt.subjs_all{aa});
    savename = fullfile(saveDir,name);
    save(savename,'cnt');
end


%% Simulates noise data based on specified alpha value

mrk = loadData(opt.subjs_all{1},'NoisePer');
newtimes = convTimes(mrk.time,100);

excludeFinal = 1000;

% gets standard deviation of channel Cz of EEG data (average std of phase 1
% and phase RT)
stds = nan(Ns,2);
for aa = 1:Ns
    
    % phase 1
    data = [];
    
    % loads data
    [mrk,cnt] = loadData(opt.subjs_all{aa},'Phase1');
    cnt = proc_selectChannels(cnt,'Cz'); % uses only data from Cz
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset',...
        'trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,2000,{'trial start','movement onset'});
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset'});
    
    % isolates within-trial data (excluding the RP period)
    times = convTimes(mrk.time,100);
    idx = ones(size(times));
    idx(1:2:end) = 0;
    times(logical(idx)) = times(logical(idx)) - (excludeFinal/10);
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    stds(aa,1) = std(data); % gets standard deviation of data
    
    % phase RT
    data = [];
    
    % loads data
    [mrk,cnt] = loadData(opt.subjs_all{aa},'RT');
    cnt = proc_selectChannels(cnt,'Cz'); % uses only data from Cz
    mrk = mrk_selectClasses(mrk,{'trial start','go signal',...
        'trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,1000,{'trial start','go signal'});
    mrk = mrk_selectClasses(mrk,{'trial start','go signal'});
    
    % isolates within-trial data (excluding the RP period)
    times = convTimes(mrk.time,100);
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    stds(aa,2) = std(data); % gets standard deviation of data
    
end
stdVals = mean(stds,2);

% sets alpha value of noise to be simulated
if alpha == 0
    app = '000';
elseif alpha == 0.5
    app = '050';
elseif alpha == 1
    app = '100';
elseif alpha == 1.5
    app = '150';
elseif alpha == 1.75
    app = '175';
else
    error('unrecognised alpha value')
end


saveDir = ''; % filepath to save simulated noise
for aa = 1:Ns
    % generates noise
    base = single(randn(newtimes(end),1));
    noiseOut = genNoise(alpha,base);
    % scales noise
    scalednoise = noiseOut/std(noiseOut);
    scalednoise = scalednoise*stdVals(aa);
    % puts into cnt format and saves cnt
    cnt.x = scalednoise;
    cnt.clab = {'Cz'};
    name = sprintf('cnt_%s_Cz_%s',opt.subjs_all{aa},app);
    savename = fullfile(saveDir,name);
    save(savename,'cnt');
end


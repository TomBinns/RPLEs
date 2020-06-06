%% Tests the prediction that low-frequencies are more dominant in voluntary
% action conditions compared to instructed action conditions

% Author(s): Thomas Binns, 2020

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

channels = {'FC2','FC1','C1','Cz','C2','CP1','CPz','CP2'}; % central-most
    % electrodes used to test prediction
         
%% gets the PSDs for phase 1 and phase RT EEG data of each participant

excludeFinal = 1000; % excludes RP period from phase 1 data
Hz = 30; % frequencies up to this value will be analysed
alphas = nan(Ns,2);
fits = cell(Ns,2);
yVals = nan(Ns,Hz,2);
h = cell(1,Ns);
for aa = 1:Ns
    
    % phase 1
    data = [];
    
    [mrk,cnt] = loadData(opt.subjs_all{aa},'Phase1');
    cnt = proc_selectChannels(cnt,channels);
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset',...
        'trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,2000,{'trial start','movement onset'});
    mrk = mrk_selectClasses(mrk,{'trial start','movement onset'});
    
    % isolates within-trial data
    times = convTimes(mrk.time,100);
    idx = ones(size(times));
    idx(1:2:end) = 0;
    times(logical(idx)) = times(logical(idx)) - (excludeFinal/10);
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    cnt.x = data;
    
    % gets alpha value of EEG data
    spec = proc_spectrum(cnt,[1,Hz],'Scaling','power');
    spec.x = mean(spec.x,2);
    xvals = log10(spec.t);
    yvals = log10(spec.x);
    yvals = yvals';
    fit = polyfit(xvals,yvals,1);
    alphas(aa,1) = fit(1)*-1;
    fits{aa,1} = fit;
    yVals(aa,:,1) = yvals;
    
    % phase RT
    data = [];
    
    [mrk,cnt] = loadData(opt.subjs_all{aa},'RT');
    cnt = proc_selectChannels(cnt,channels);
    mrk = mrk_selectClasses(mrk,{'trial start','go signal','trial end'});
    trial_mrk = getTrialMarkers(mrk);
    trial_mrk = trial_mrk(cellfun(@length,trial_mrk)==3);
    mrk = mrk_selectEvents(mrk,[trial_mrk{:}]);
    mrk = removeShortTrials(mrk,1000,{'trial start','go signal'});
    mrk = mrk_selectClasses(mrk,{'trial start','go signal'});
    
    % isolates data from trials
    times = convTimes(mrk.time,100);
    idx = ones(size(times));
    idx(1:2:end) = 0;
    for bb = 1:length(times)/2
        data = cat(1,data,cnt.x(times((bb*2)-1):times(bb*2),:));
    end
    cnt.x = data;
    
    % gets alpha value of subject's EEG data
    spec = proc_spectrum(cnt,[1,Hz],'Scaling','power');
    spec.x = mean(spec.x,2);
    xvals = log10(spec.t);
    yvals = log10(spec.x);
    yvals = yvals';
    fit = polyfit(xvals,yvals,1);
    alphas(aa,2) = fit(1)*-1;
    fits{aa,2} = fit;
    yVals(aa,:,2) = yvals;
    
    
    % plots PSDs with alpha values
    if aa == 1 || aa == 10
        figure
        hold on
    end
    if aa < 10
        subplot(3,3,aa)
        hold on
        grid on
    end
    if aa > 9
        subplot(3,3,aa-9)
        hold on
        grid on
    end
    h{aa} = plot((1:Hz),yVals(aa,:,1),'Color',[0,0,0],'LineWidth',3);
    plot((1:Hz),yVals(aa,:,2),'Color',[.4,.4,.4],'LineWidth',3);
    ttl = sprintf('Phase 1 \\alpha = %.2f\nPhase RT \\alpha =  %.2f',...
        alphas(aa,:));
    title(ttl);
    txt = sprintf('Participant %.0f',aa);
    xCoord = (Hz/100)*65;
    yCoord = h{aa}.Parent.YLim(1)+...
        (((h{aa}.Parent.YLim(2)-h{aa}.Parent.YLim(1))/100)*90);
    text(xCoord,yCoord,txt,'FontSize',15,'FontWeight','bold');
    
    h{aa}.Parent.FontSize = 15;
    h{aa}.Parent.FontWeight = 'bold';
    h{aa}.Parent.XColor = [0,0,0];
    h{aa}.Parent.YColor = [0,0,0];
    
    if aa == 1
        legend('Phase 1','Phase RT');
    end    
    if aa == 10
        ylabel('Log_{10} Power (\muV^2/Hz)');
        h{aa}.Parent.YLabel.FontSize = 20;
    end
    if aa == Ns
        xlabel('Frequency (Hz)');
        h{aa}.Parent.XLabel.FontSize = 20;
    end
    
end

meanalphas = mean(alphas,2);


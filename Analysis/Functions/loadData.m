function [mrk,cnt,mnt] = loadData(subj_code,phase_name)
% loadData - Loads EEG or noise data, trial markers, and electrode data
%     (latter for EEG data only)
% 
% Argument(s):
%  subj_code -    identifying code of subject whose data is to be loaded 
%                 (e.g. 'VPtab')
%  phase_name -   specifies which experimental phase should be loaded; 
%                 'Phase1' for phase 1, 'RT' for phase RT, 'NoiseXXX' for
%                 simulated noise, where XXX specifies the alpha value of
%                 the noise (e.g. 000 = 0.00, 150 = 1.50, Per = personal)
% 
% Returns:
%  mrk -          trial markers
%  cnt -          recorded/simulated EEG data
%  mnt -          electrodes & their positions (for EEG data only)
%  
% Author(s): Matthias Schultze-Kraft; Thomas Binns, 2020


global opt

if strcmp(phase_name,'Phase1') || strcmp(phase_name,'RT')
    
    session_name = 'TrafficLight';
    
    folder = 'cnt+mrk';
    filename_eeg = sprintf('%s_%s_%s',session_name,phase_name,subj_code);
    filename_eeg = fullfile(opt.dataDir,folder,filename_eeg);
    filename_mrk = sprintf('%s_mrk.mat',filename_eeg);
    
    fprintf('Loading data set %s, %s...\n',subj_code,phase_name)
    
    [cnt,mrk,mnt] = file_loadMatlab(filename_eeg);
    if exist(filename_mrk,'file')
        load(filename_mrk)
    end
    
    cnt = proc_selectChannels(cnt,opt.clab_load);
    mnt = mnt_restrictMontage(mnt,cnt);
    mnt.scale_box = [];
    mnt = mnt_scalpToGrid(mnt);
    
    if strcmp(phase_name,'Phase1')
        switch session_name
            case 'TrafficLight'
                ci = strcmp(mrk.className,'start phase1');
                mrk.className{ci} = 'trial start';
                ci = strcmp(mrk.className,'EMG onset');
                mrk.className{ci} = 'movement onset';
        end
    elseif strcmp(phase_name,'RT')
        switch session_name
            case 'TrafficLight'
                ci = strcmp(mrk.className,'start rt');
                mrk.className{ci} = 'trial start';
                ci = strcmp(mrk.className,'EMG onset');
                mrk.className{ci} = 'movement onset';
                ci = strcmp(mrk.className,'light rt');
                mrk.className{ci} = 'go signal';
        end
    end

elseif strcmp(phase_name(1:5),'Noise')
    
    folder = 'cnt+mrk';
    filename_eeg = sprintf('cnt_%s_Cz_%s.mat',subj_code,phase_name(6:8));
    filename_eeg = fullfile(opt.dataDir,folder,filename_eeg);
    load(filename_eeg);
    filename_mrk = fullfile(opt.dataDir,'cnt+mrk\noisemrk.mat');
    load(filename_mrk);
    mnt = [];
    
else
    
    error('Unknown phase requested')
    
end


end
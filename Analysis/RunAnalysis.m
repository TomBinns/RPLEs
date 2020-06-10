%% Sets up analysis

global opt

% specify directory where data is located
opt.dataDir = 'C:\Users\tomth\Documents\Coding\RPLE_data';

%% Gets and plots results for EEG data (number of events, similarity 
% scores, appearance of events); may take a while!

[stats_NRPLEs,stats_SimScores] = analysis_EEG;

%% Gets and plots results for noise data (number of events, similarity 
% scores, appearance of events); may take a while!

[stats_Noise] = analysis_noise;


function trial_mrk = getTrialMarkers(mrk,must_contain)
% getTrialMarkers - shows which markers correspond to each trial
% 
% Argument(s):
%  mrk -      trial markers
% 
% Returns:
%  trial_mrk - cell for each trial containing the markers which correspond 
%              to that trial, numbered based on which column they are found 
%              in mrk.y
%              
% Author(s): Matthias-Schultze Kraft             

mrk = mrk_sortChronologically(mrk);
ci_ts = find(strcmp(mrk.className,'trial start'));
ci_te = find(strcmp(mrk.className,'trial end'));

% extract markers
idx = 1;
trial_mrk = [];
n_trial = 0;
while idx<length(mrk.time)
    if find(mrk.y(:,idx))==ci_ts
        event_idx = [];
        class_idx = [];
        while 1
            event_idx = [event_idx idx];
            class_idx = [class_idx find(mrk.y(:,idx))];
            if find(mrk.y(:,idx))==ci_te
                break
            end
            idx = idx+1;
        end
        n_trial = n_trial+1;
        trial_mrk{n_trial} = event_idx;
    end
    idx = idx+1;
end

if nargin>1
    if ischar(must_contain)
        must_contain = {must_contain};
    end
    trial_mrk = select_trials(trial_mrk,mrk,must_contain);
end


function trial = select_trials(trial,mrk,must_contain)

Nc = length(must_contain);
not_mode = zeros(Nc,1);
for kk = 1:Nc
    if strcmp(must_contain{kk}(1),'~')
        must_contain{kk}(1) = [];
        not_mode(kk) = 1;
    end    
end
Nt = length(trial);
keep = [];
for jj = 1:Nt
    mrk2 = mrk_selectEvents(mrk,trial{jj});
    conform = zeros(Nc,1);
    for kk = 1:Nc
        conform(kk) = ismember(must_contain{kk},mrk2.className);
        if not_mode(kk)
            conform(kk) = ~conform(kk);
        end
    end
    if all(conform)
        keep = [keep jj];
    end
end
trial = trial(keep);


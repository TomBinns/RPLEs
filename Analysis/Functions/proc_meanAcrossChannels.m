function epo= proc_meanAcrossChannels(epo)
%PROC_MEANACROSSCHANNELS - Average signals across all channels
%
%Synopsis:
%  EPO= proc_meanAcrossChannels(EPO)
%
%Arguments:
%  DAT  - data structure of epoched data
%
%Returns:
%  DAT  - updated data structure


epo.x = mean(epo.x(:,:,:),2);
epo.clab = {'1'};


end


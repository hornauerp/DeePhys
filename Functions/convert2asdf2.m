function spks_data = convert2asdf2(spike_times, cluster_ids, sampling_rate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spks_data.binsize = 1/(sampling_rate/1000);
spks_data.nbins = round(sampling_rate*max(spike_times));
spks_data.nchannels = length(unique(cluster_ids));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ch_unique = unique(cluster_ids);
for j = 1:spks_data.nchannels
    spk_time_j = spike_times(cluster_ids == ch_unique(j)); 
    spks_data.raster(j,:) = {round(spk_time_j*sampling_rate)'};
end
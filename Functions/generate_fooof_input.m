function data = generate_fooof_input(bursts,units,ops)
arguments
    bursts cell
    units cell %if left empty, will be calculated on the whole network instead of individual units
    ops.binning double = 0.001 %in seconds
    ops.min_length double = 0 %in seconds
    ops.min_participation double = 0.2 %minimum ratio of bursts a unit is in to be considered for analysis
    ops.concatenated logical = false %Flag whether bursts should be concatenated or not
    % only relevant if concatenated=true
    ops.duration double = 60 %we split up the concatenated signal into chunks of this length
    ops.overlap double = 0.1 %with this much overlap
end

if isempty(units)
    units = cellfun(@(x) ones(size(x)),bursts,'un',0);
    good_unit_ids = 1;
else
    unit_ids = cellfun(@(x) unique(x)',units,'un',0);
    [g,gc] = histcounts([unit_ids{:}],1:max([unit_ids{:}]));

    good_unit_ids = gc(g>length(units)*ops.min_participation);
end


if ops.concatenated
    long_b = vertcat(bursts{:});
    long_b = {smoothdata(histcounts(long_b,long_b(1):ops.binning:long_b(end)),"gaussian",3)};
    total_duration = size(long_b{1},2) * ops.binning;
    % total_duration = size(long_b{1},2) * ops.binning;
    % fprintf('Total burst duration: %.3f s\n',total_duration)
else
    burst_matrices = cell(1,length(bursts));
    for iburst = 1:length(bursts)
        burst_duration = round((bursts{iburst}(end) - bursts{iburst}(1)) / ops.binning);
        burst_mat = zeros(length(good_unit_ids), burst_duration);
        for iunit = 1:length(good_unit_ids)
            unit_spikes = bursts{iburst}(units{iburst}==good_unit_ids(iunit));
            % unit_spikes = unit_spikes - unit_spikes(1);
            % if length(unit_spikes) > 5
            ft_input = histcounts(unit_spikes,unit_spikes(1):ops.binning:unit_spikes(1)+(burst_duration*ops.binning));
            burst_mat(iunit,:) = smoothdata(ft_input,"gaussian",3);
            % end
        end
        burst_matrices{iburst} = burst_mat;
    end

    long_b = burst_matrices(cellfun(@(x) (size(x,2)*ops.binning) > ops.min_length,burst_matrices));
end

data = struct();
for iburst = 1:length(long_b)
    ft_input = long_b{iburst};
    time = 0:ops.binning:ops.binning*size(ft_input,2); time = time(2:end);
    data.trial{1,iburst} = ft_input;
    data.time{1,iburst} = time;
end
data.label     = arrayfun(@num2str, good_unit_ids,'un',0)';
if ops.concatenated
    cfg               = [];
    cfg.length        = min([ops.duration, total_duration]) ;
    cfg.overlap       = ops.overlap;
    data              = ft_redefinetrial(cfg, data);
end
end

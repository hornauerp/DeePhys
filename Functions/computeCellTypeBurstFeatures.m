function [network_table, unit_table] = computeCellTypeBurstFeatures( ...
    bursts, spike_times, spike_units, units, cell_type_labels, unit_feature_table)
% COMPUTECELLTYPEBURSTFEATURES  Burst dynamics stratified by cell type.
%
%   [network_table, unit_table] = computeCellTypeBurstFeatures( ...
%       bursts, spike_times, spike_units, units, cell_type_labels, unit_feature_table)
%
% INPUTS:
%   bursts             - struct with .T_start and .T_end (1 x N_bursts, seconds)
%   spike_times        - (N_spikes x 1) spike times in seconds
%   spike_units        - (N_spikes x 1) template ID per spike
%   units              - (1 x N_units) UnitData array (for .TemplateID)
%   cell_type_labels   - (1 x N_units) 1=excitatory, 2=inhibitory, NaN=unclassified
%   unit_feature_table - table with BurstFractionSpikes column (optional)
%
% OUTPUTS:
%   network_table - (1 x 9) table of network-level burst features by cell type
%   unit_table    - (N_units x 1) table with BurstParticipationRate per unit

    n_units  = numel(cell_type_labels);
    n_bursts = numel(bursts.T_start);

    template_ids    = arrayfun(@(u) u.TemplateID, units);
    exc_idx         = find(cell_type_labels == 1);
    inh_idx         = find(cell_type_labels == 2);
    exc_tids        = template_ids(exc_idx);
    inh_tids        = template_ids(inh_idx);
    exc_mask_spikes = ismember(spike_units, exc_tids);
    inh_mask_spikes = ismember(spike_units, inh_tids);

    % Per-unit burst participation rate
    burst_participation = NaN(n_units, 1);
    for u = 1:n_units
        u_spikes = spike_times(spike_units == template_ids(u));
        if n_bursts == 0
            burst_participation(u) = 0;
            continue
        end
        n_participated = 0;
        for b = 1:n_bursts
            if any(u_spikes >= bursts.T_start(b) & u_spikes <= bursts.T_end(b))
                n_participated = n_participated + 1;
            end
        end
        burst_participation(u) = n_participated / n_bursts;
    end
    unit_table = table(burst_participation, 'VariableNames', {'BurstParticipationRate'});

    if n_bursts == 0
        nw = makeNaNStruct();
        network_table = struct2table(nw);
        return
    end

    % Per-burst first-spike latencies and type-level participation
    lead_times_E = NaN(1, n_bursts);
    lead_times_I = NaN(1, n_bursts);
    has_E        = false(1, n_bursts);
    has_I        = false(1, n_bursts);
    e_led        = false(1, n_bursts);

    for b = 1:n_bursts
        t0   = bursts.T_start(b);
        t1   = bursts.T_end(b);
        in_b = spike_times >= t0 & spike_times <= t1;

        e_t = spike_times(in_b & exc_mask_spikes);
        i_t = spike_times(in_b & inh_mask_spikes);

        if ~isempty(e_t)
            lead_times_E(b) = min(e_t) - t0;
            has_E(b)        = true;
        end
        if ~isempty(i_t)
            lead_times_I(b) = min(i_t) - t0;
            has_I(b)        = true;
        end
        if has_E(b) && has_I(b)
            e_led(b) = lead_times_E(b) < lead_times_I(b);
        end
    end

    both = has_E & has_I;
    BurstLeadTime_E   = mean(lead_times_E, 'omitnan');
    BurstLeadTime_I   = mean(lead_times_I, 'omitnan');
    if sum(both) > 0
        BurstLeadFraction = sum(e_led & both) / sum(both);
    else
        BurstLeadFraction = NaN;
    end

    BurstParticipation_E = mean(has_E);
    BurstParticipation_I = mean(has_I);

    % Mean BurstFractionSpikes by type (from existing per-unit column)
    MeanBurstFraction_E = NaN;
    MeanBurstFraction_I = NaN;
    if ~isempty(unit_feature_table) && ...
            ismember('BurstFractionSpikes', unit_feature_table.Properties.VariableNames)
        bfs = unit_feature_table.BurstFractionSpikes;
        if ~isempty(exc_idx)
            MeanBurstFraction_E = mean(bfs(exc_idx), 'omitnan');
        end
        if ~isempty(inh_idx)
            MeanBurstFraction_I = mean(bfs(inh_idx), 'omitnan');
        end
    end

    % Intra-burst firing rates
    total_burst_dur = sum(bursts.T_end - bursts.T_start);
    IntraBurstRate_E = NaN;
    IntraBurstRate_I = NaN;
    if total_burst_dur > 0
        in_any_burst = false(size(spike_times));
        for b = 1:n_bursts
            in_any_burst = in_any_burst | ...
                (spike_times >= bursts.T_start(b) & spike_times <= bursts.T_end(b));
        end
        n_exc = numel(exc_idx);
        n_inh = numel(inh_idx);
        if n_exc > 0
            IntraBurstRate_E = sum(in_any_burst & exc_mask_spikes) / (n_exc * total_burst_dur);
        end
        if n_inh > 0
            IntraBurstRate_I = sum(in_any_burst & inh_mask_spikes) / (n_inh * total_burst_dur);
        end
    end

    nw = struct( ...
        'BurstLeadTime_E',      BurstLeadTime_E, ...
        'BurstLeadTime_I',      BurstLeadTime_I, ...
        'BurstLeadFraction',    BurstLeadFraction, ...
        'MeanBurstFraction_E',  MeanBurstFraction_E, ...
        'MeanBurstFraction_I',  MeanBurstFraction_I, ...
        'IntraBurstRate_E',     IntraBurstRate_E, ...
        'IntraBurstRate_I',     IntraBurstRate_I, ...
        'BurstParticipation_E', BurstParticipation_E, ...
        'BurstParticipation_I', BurstParticipation_I);
    network_table = struct2table(nw);
end


function nw = makeNaNStruct()
    nw = struct( ...
        'BurstLeadTime_E',      NaN, ...
        'BurstLeadTime_I',      NaN, ...
        'BurstLeadFraction',    NaN, ...
        'MeanBurstFraction_E',  NaN, ...
        'MeanBurstFraction_I',  NaN, ...
        'IntraBurstRate_E',     NaN, ...
        'IntraBurstRate_I',     NaN, ...
        'BurstParticipation_E', NaN, ...
        'BurstParticipation_I', NaN);
end

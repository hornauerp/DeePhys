function [network_table, unit_table] = computeSpatialBurstFeatures( ...
    electrode_coordinates, units, bursts, spike_times, spike_units, cell_type_labels, params)
% COMPUTESPATIALBURSTFEATURES  Spatial origin and propagation of network bursts.
%
%   [network_table, unit_table] = computeSpatialBurstFeatures( ...
%       electrode_coordinates, units, bursts, spike_times, spike_units, cell_type_labels, params)
%
% For each network burst, identifies the spatial origin (first-spike location)
% and estimates the propagation wavefront speed via regression of first-spike
% latency on distance from the burst origin.
%
% INPUTS:
%   electrode_coordinates - (N_channels x 2) XY positions in micrometers
%   units                 - (1 x N_units) UnitData array with .ReferenceElectrode, .TemplateID
%   bursts                - struct with .T_start and .T_end (1 x N_bursts, seconds)
%   spike_times           - (N_spikes x 1) spike times in seconds
%   spike_units           - (N_spikes x 1) template ID per spike
%   cell_type_labels      - (1 x N_units) 1=exc, 2=inh, NaN=unclassified (may be [])
%   params                - struct with optional field:
%                             .BurstMinUnitsForWavefront (default 5)
%
% OUTPUTS:
%   network_table - (1 x 3) table: BurstOriginDispersion, MeanBurstPropagationSpeed,
%                                  BurstOriginExcFraction
%   unit_table    - empty table (no per-unit features added here)

    arguments
        electrode_coordinates (:,2) double
        units                 (1,:)
        bursts                struct
        spike_times           (:,1) double
        spike_units           (:,1) double
        cell_type_labels      (1,:) double = []
        params                struct = struct()
    end

    p = mergeDefaults(params, struct('BurstMinUnitsForWavefront', 5));

    nw.BurstOriginDispersion      = NaN;
    nw.MeanBurstPropagationSpeed  = NaN;
    nw.BurstOriginExcFraction     = NaN;
    network_table = struct2table(nw);
    unit_table    = table();

    n_units  = numel(units);
    n_bursts = numel(bursts.T_start);

    if n_bursts == 0 || n_units < p.BurstMinUnitsForWavefront
        return
    end

    ref_els      = arrayfun(@(u) u.ReferenceElectrode, units);
    template_ids = arrayfun(@(u) u.TemplateID, units);
    xy           = electrode_coordinates(ref_els, :);   % (N_units x 2)

    has_labels = ~isempty(cell_type_labels) && numel(cell_type_labels) == n_units;

    % Precompute per-unit spike times to avoid O(n_bursts x n_units x n_spikes) scan
    unit_spikes = cell(n_units, 1);
    for u = 1:n_units
        unit_spikes{u} = spike_times(spike_units == template_ids(u));
    end

    % Per-burst analysis
    origin_xy     = NaN(n_bursts, 2);
    prop_speeds   = NaN(n_bursts, 1);
    origin_is_exc = NaN(n_bursts, 1);

    for b = 1:n_bursts
        t0 = bursts.T_start(b);
        t1 = bursts.T_end(b);

        % First-spike time per unit during this burst
        first_t = NaN(n_units, 1);
        for u = 1:n_units
            st_u   = unit_spikes{u};
            u_in_b = st_u >= t0 & st_u <= t1;
            if any(u_in_b)
                first_t(u) = min(st_u(u_in_b)) - t0;  % latency from burst start
            end
        end

        participating = find(~isnan(first_t));
        if numel(participating) < p.BurstMinUnitsForWavefront, continue; end

        % Burst origin: unit with earliest first spike
        [~, origin_unit] = min(first_t);
        origin_xy(b,:)   = xy(origin_unit, :);

        if has_labels && ~isnan(cell_type_labels(origin_unit))
            origin_is_exc(b) = cell_type_labels(origin_unit) == 1;
        end

        % Propagation speed: regress latency ~ distance from origin
        dists = sqrt(sum((xy(participating,:) - origin_xy(b,:)).^2, 2));
        lats  = first_t(participating);

        % Only use units with positive distance from origin
        moving = dists > 0;
        if sum(moving) < 3, continue; end

        % slope: dt/dx (s/um), speed = 1/slope (um/s -> um/ms * 0.001)
        p_fit = polyfit(dists(moving), lats(moving) * 1000, 1);  % latency in ms
        if p_fit(1) > 0
            prop_speeds(b) = 1 / p_fit(1);  % um/ms
        end
    end

    valid_origins = ~any(isnan(origin_xy), 2);
    if sum(valid_origins) >= 2
        centroid = mean(origin_xy(valid_origins,:), 1);
        nw.BurstOriginDispersion = mean( ...
            sqrt(sum((origin_xy(valid_origins,:) - centroid).^2, 2)), 'omitnan');
    end

    if any(~isnan(prop_speeds))
        nw.MeanBurstPropagationSpeed = mean(prop_speeds, 'omitnan');
    end

    valid_type = ~isnan(origin_is_exc);
    if sum(valid_type) > 0
        nw.BurstOriginExcFraction = sum(origin_is_exc(valid_type)) / sum(valid_type);
    end

    network_table = struct2table(nw);
end


function p = mergeDefaults(p, defaults)
    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        if ~isfield(p, fields{i})
            p.(fields{i}) = defaults.(fields{i});
        end
    end
end

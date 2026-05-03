function ax = plotSpatialUnitMap(electrode_coordinates, units, cell_type_labels, options)
% PLOTSPATIALUNITMAP  Scatter map of units at their reference electrode positions.
%
%   ax = plotSpatialUnitMap(electrode_coordinates, units, cell_type_labels)
%   ax = plotSpatialUnitMap(..., Name, Value)
%
% Plots each unit as a filled circle at its reference electrode XY position,
% colored by cell type. Optionally scales marker size by firing rate.
%
% INPUTS:
%   electrode_coordinates - (N_channels x 2) XY positions in micrometers
%   units                 - (1 x N_units) UnitData array with .ReferenceElectrode
%   cell_type_labels      - (1 x N_units) 1=excitatory, 2=inhibitory, NaN=unclassified
%                           (pass [] to skip cell-type coloring)
%
% NAME-VALUE OPTIONS:
%   FiringRates  - (N_units x 1) firing rates for marker size scaling (optional)
%   ShowGrid     - logical (default true) show all electrode positions as background
%   Ax           - target axes (creates new figure if empty)
%   MarkerSize   - base marker size (default 60)
%   Title        - axes title string
%
% OUTPUT:
%   ax - axes handle

    arguments
        electrode_coordinates (:,2) double
        units                 (1,:)
        cell_type_labels      (1,:) double = []
        options.FiringRates   (:,1) double = []
        options.ShowGrid      (1,1) logical = true
        options.Ax                         = []
        options.MarkerSize    (1,1) double = 60
        options.Title         (1,:) char   = ''
    end

    % Colors: excitatory=blue, inhibitory=red, unclassified=grey
    c_exc  = [0.22, 0.49, 0.72];
    c_inh  = [0.80, 0.17, 0.15];
    c_unk  = [0.65, 0.65, 0.65];

    if isempty(options.Ax)
        figure('Color', 'w');
        ax = axes;
    else
        ax = options.Ax;
    end
    hold(ax, 'on');

    % Background electrode grid
    if options.ShowGrid
        scatter(ax, electrode_coordinates(:,1), electrode_coordinates(:,2), ...
            8, 'filled', 'Marker', 'square', ...
            'MarkerFaceAlpha', 0.12, 'MarkerFaceColor', [0.5 0.5 0.5], ...
            'MarkerEdgeAlpha', 0);
    end

    n_units = numel(units);
    if n_units == 0
        axis(ax, 'equal'); axis(ax, 'off');
        return
    end

    ref_els = arrayfun(@(u) u.ReferenceElectrode, units);
    xy      = electrode_coordinates(ref_els, :);

    % Assign colors
    colors = repmat(c_unk, n_units, 1);
    if ~isempty(cell_type_labels) && numel(cell_type_labels) == n_units
        colors(cell_type_labels == 1, :) = repmat(c_exc, sum(cell_type_labels == 1), 1);
        colors(cell_type_labels == 2, :) = repmat(c_inh, sum(cell_type_labels == 2), 1);
    end

    % Marker sizes (scaled by firing rate if provided)
    base_sz = options.MarkerSize;
    if ~isempty(options.FiringRates) && numel(options.FiringRates) == n_units
        fr = options.FiringRates(:);
        fr_norm = (fr - min(fr)) / max(max(fr) - min(fr), eps);
        sizes = base_sz * (0.3 + 1.5 * fr_norm);
    else
        sizes = base_sz * ones(n_units, 1);
    end

    % Plot each unit
    scatter(ax, xy(:,1), xy(:,2), sizes, colors, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');

    % Legend entries
    h = gobjects(0);
    lbl = {};
    if ~isempty(cell_type_labels) && numel(cell_type_labels) == n_units
        n_exc = sum(cell_type_labels == 1);
        n_inh = sum(cell_type_labels == 2);
        n_unk = sum(isnan(cell_type_labels));
        if n_exc > 0
            h(end+1) = scatter(ax, NaN, NaN, base_sz, c_exc, 'filled');
            lbl{end+1} = sprintf('Excitatory (n=%d)', n_exc);
        end
        if n_inh > 0
            h(end+1) = scatter(ax, NaN, NaN, base_sz, c_inh, 'filled');
            lbl{end+1} = sprintf('Inhibitory (n=%d)', n_inh);
        end
        if n_unk > 0
            h(end+1) = scatter(ax, NaN, NaN, base_sz, c_unk, 'filled');
            lbl{end+1} = sprintf('Unclassified (n=%d)', n_unk);
        end
        legend(ax, h, lbl, 'Location', 'best', 'Box', 'off');
    end

    xlabel(ax, 'x (\mum)');
    ylabel(ax, 'y (\mum)');
    if ~isempty(options.Title)
        title(ax, options.Title);
    end
    axis(ax, 'equal');
    ax.Box = 'off';
    hold(ax, 'off');
end

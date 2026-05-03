function ax = plotSpatialEIBalance(viz_data, electrode_coordinates, options)
% PLOTSPATIALEIBALANCE  Heatmap of local excitatory fraction on the MEA chip.
%
%   ax = plotSpatialEIBalance(viz_data, electrode_coordinates)
%   ax = plotSpatialEIBalance(viz_data, electrode_coordinates, Name, Value)
%
% Visualizes the output of computeSpatialEIBalance as a heatmap overlaid
% on the electrode grid. Color encodes local excitatory fraction (0=inhibitory,
% 1=excitatory, 0.5=balanced, shown as diverging colormap).
%
% INPUTS:
%   viz_data              - struct returned by computeSpatialEIBalance with fields:
%                             .grid_x, .grid_y, .exc_fraction, .unit_xy, .unit_labels
%   electrode_coordinates - (N_channels x 2) XY positions in micrometers
%
% NAME-VALUE OPTIONS:
%   ShowUnits    - logical (default true) overlay unit positions as dots
%   ShowGrid     - logical (default false) show electrode grid
%   Ax           - target axes
%   Title        - axes title string

    arguments
        viz_data              struct
        electrode_coordinates (:,2) double = zeros(0,2)
        options.ShowUnits     (1,1) logical = true
        options.ShowGrid      (1,1) logical = false
        options.Ax                         = []
        options.Title         (1,:) char   = 'Local E/I balance'
    end

    required = {'grid_x', 'grid_y', 'exc_fraction'};
    for r = required
        if ~isfield(viz_data, r{1})
            error('plotSpatialEIBalance: viz_data missing field "%s"', r{1});
        end
    end

    if isempty(options.Ax)
        figure('Color', 'w');
        ax = axes;
    else
        ax = options.Ax;
    end
    hold(ax, 'on');

    % Heatmap
    pcolor(ax, viz_data.grid_x, viz_data.grid_y, viz_data.exc_fraction);
    shading(ax, 'interp');
    colormap(ax, buildDivergingColormap());
    cb = colorbar(ax);
    cb.Label.String = 'Local excitatory fraction';
    caxis(ax, [0 1]);

    if options.ShowGrid && ~isempty(electrode_coordinates)
        scatter(ax, electrode_coordinates(:,1), electrode_coordinates(:,2), ...
            6, 'filled', 'Marker', 'square', ...
            'MarkerFaceAlpha', 0.15, 'MarkerFaceColor', [0.3 0.3 0.3], ...
            'MarkerEdgeAlpha', 0);
    end

    if options.ShowUnits && isfield(viz_data, 'unit_xy') && ~isempty(viz_data.unit_xy)
        xy = viz_data.unit_xy;
        c_exc = [0.22, 0.49, 0.72];
        c_inh = [0.80, 0.17, 0.15];
        c_unk = [0.65, 0.65, 0.65];

        if isfield(viz_data, 'unit_labels')
            lbl = viz_data.unit_labels(:);
            colors = repmat(c_unk, size(xy,1), 1);
            colors(lbl == 1, :) = repmat(c_exc, sum(lbl == 1), 1);
            colors(lbl == 2, :) = repmat(c_inh, sum(lbl == 2), 1);
        else
            colors = repmat(c_unk, size(xy,1), 1);
        end
        scatter(ax, xy(:,1), xy(:,2), 25, colors, 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
    end

    xlabel(ax, 'x (\mum)');
    ylabel(ax, 'y (\mum)');
    title(ax, options.Title);
    axis(ax, 'equal');
    ax.Box = 'off';
    hold(ax, 'off');
end


function cmap = buildDivergingColormap()
% Red (inhibitory) - white (balanced) - blue (excitatory) diverging colormap.
    n    = 64;
    half = n / 2;
    r_to_w = [linspace(0.80, 1, half)', linspace(0.17, 1, half)', linspace(0.15, 1, half)'];
    w_to_b = [linspace(1, 0.22, half)', linspace(1, 0.49, half)', linspace(1, 0.72, half)'];
    cmap   = [r_to_w; w_to_b];
end

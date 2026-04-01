function ax = PlotNetworkActivity(eia, culture_idx, opts)
% PLOTNETWORKACTIVITY  Plot total network firing rate coloured by E/I burst state.
%
% The firing rate trace is coloured per bin: burst-state bins (excitatory-dominant)
% are red, non-burst bins are blue.
%
% INPUTS:
%   eia          - EIAnalyzer object
%   culture_idx  - Index of the culture to plot (default 1)
%
% OPTIONS (name-value):
%   Parent       - axes to plot into; creates a new figure if not provided
%
% OUTPUTS:
%   ax           - axes handle

arguments
    eia         EIAnalyzer
    culture_idx (1,1) double = 1
    opts.Parent = []
end

assert(~isempty(eia.BurstState), ...
    'Run computeActivity() and detectBursts() before plotting');

act = eia.Activity(culture_idx);
z   = eia.BurstState{culture_idx};

if isempty(opts.Parent)
    figure('Color', 'w');
    ax = axes();
else
    ax = opts.Parent;
end

colors = [0.12, 0.47, 0.71; 0.84, 0.19, 0.15];  % blue (non-burst) / red (burst)
vertex_colors = colors(z, :);

patch(ax, ...
    'XData',           [act.x, NaN], ...
    'YData',           [act.total, NaN], ...
    'ZData',           [zeros(size(act.x)), NaN], ...
    'FaceVertexCData', [vertex_colors; [NaN NaN NaN]], ...
    'FaceColor',       'none', ...
    'EdgeColor',       'interp', ...
    'LineWidth',       0.5);
colormap(ax, colors);
xlabel(ax, 'Time (s)');
ylabel(ax, 'Network firing rate (spikes / bin)');
title(ax, sprintf('Culture %i — network activity (burst state coloured)', culture_idx), ...
    'FontWeight', 'normal');
box(ax, 'off');

setAllFontSizes(ancestor(ax, 'figure'), eia.Parameters.Plotting.FontSize);
end

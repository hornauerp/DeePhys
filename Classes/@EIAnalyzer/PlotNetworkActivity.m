function PlotNetworkActivity(eia, culture_idx)
% PLOTNETWORKACTIVITY  Plot total network firing rate coloured by E/I burst state.
%
% The firing rate trace is coloured per bin: burst-state bins (excitatory-dominant)
% are coloured differently from non-burst bins, using an RdBu colormap.
%
% INPUTS:
%   eia          - EIAnalyzer object
%   culture_idx  - Index of the culture to plot (default 1)

arguments
    eia         EIAnalyzer
    culture_idx (1,1) double = 1
end

assert(~isempty(eia.BurstState), ...
    'Run computeActivity() and detectBursts() before plotting');

act = eia.Activity(culture_idx);
z   = eia.BurstState{culture_idx};

colors = [0.12, 0.47, 0.71; 0.84, 0.19, 0.15];  % blue (burst) / red (non-burst)
vertex_colors = colors(z, :);

smooth_total = smoothdata(act.total, 'SmoothingFactor', 0.1);

figure('Color', 'w');
patch('XData',         [act.x, NaN], ...
      'YData',         [smooth_total, NaN], ...
      'ZData',         [zeros(size(act.x)), NaN], ...
      'FaceVertexCData', [vertex_colors; [NaN NaN NaN]], ...
      'FaceColor',     'none', ...
      'EdgeColor',     'interp', ...
      'LineWidth',     0.5);
colormap(colors);
xlabel('Time (s)');
ylabel('Network firing rate (spikes / bin)');
title(sprintf('Culture %i — network activity (burst state coloured)', culture_idx), ...
    'FontWeight', 'normal');
box off;
end

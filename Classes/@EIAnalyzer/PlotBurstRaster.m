function ax = PlotBurstRaster(eia, culture_idx, opts)
% PLOTBURSTRASTER  Raster plot of unit spike times aligned to accepted burst peaks.
%
% For each accepted burst, plots individual unit spike times relative to the burst
% peak, coloured by cell type (blue = excitatory, red = inhibitory, grey = unclassified).
% Units are sorted first by cell type (E → I → unclassified), then within each type
% by mean spike count during burst windows (descending).
%
% The window displayed is ±PeakCutout bins on each side of the burst peak
% (matching the burst cutout width used during extraction).
%
% INPUTS:
%   eia          - EIAnalyzer object
%   culture_idx  - Index of the culture to plot (default 1)
%
% OPTIONS (name-value):
%   BurstIndex   - Scalar index to plot a single burst, or 'all' to overlay all
%                  accepted bursts (default 'all')
%   ClusterIndex - Plot only bursts belonging to this cluster ([] = all clusters,
%                  default [])
%   Parent       - axes to plot into; creates a new figure if not provided
%
% OUTPUTS:
%   ax           - axes handle
%
% Requires eia.Activity and eia.BurstCutouts (run computeActivity, detectBursts,
% and extractBurstCutouts first).

arguments
    eia         EIAnalyzer
    culture_idx (1,1) double = 1
    opts.BurstIndex  = 'all'  % 'all' | scalar integer
    opts.ClusterIndex (1,:) double = []
    opts.Parent = []
end

assert(~isempty(eia.BurstCutouts), ...
    'Run extractBurstCutouts() before PlotBurstRaster()');
assert(culture_idx >= 1 && culture_idx <= numel(eia.Activity), ...
    'culture_idx %d out of range [1, %d]', culture_idx, numel(eia.Activity));

bc  = eia.BurstCutouts(culture_idx);
act = eia.Activity(culture_idx);
p   = eia.Parameters.BurstDetection;

if isempty(bc.accepted.total)
    warning('EIAnalyzer:PlotBurstRaster', ...
        'No accepted bursts for culture %d — skipping.', culture_idx);
    ax = gobjects(0);
    return
end

% ── Unit info for this culture ────────────────────────────────────────────
ctc        = eia.Classifier;
map        = eia.getCultureUnitMapping();
cid        = map.unique_cultures(culture_idx);
unit_mask  = map.culture_ids == cid;
unit_rows  = find(unit_mask);
ud_indices = map.ud_order(unit_rows);
valid      = ud_indices > 0;
unit_rows  = unit_rows(valid);
ud_indices = ud_indices(valid);
labels     = map.labels(unit_rows);   % 1=exc, 2=inh, NaN=unclassified

ud     = ctc.UnitDataArray;
bin_sz = eia.Parameters.Activity.BinSize;  % seconds per bin
half_w = p.PeakCutout * bin_sz;            % half-window in seconds

% ── Select burst peaks ────────────────────────────────────────────────────
peak_locs = bc.accepted_locs;  % bin indices of accepted burst peaks

if ~isempty(opts.ClusterIndex) && ~isempty(bc.cluster_idx)
    cluster_mask = ismember(bc.cluster_idx, opts.ClusterIndex);
    peak_locs    = peak_locs(cluster_mask);
end

if isnumeric(opts.BurstIndex) && isscalar(opts.BurstIndex)
    bi = opts.BurstIndex;
    assert(bi >= 1 && bi <= numel(peak_locs), ...
        'EIAnalyzer:PlotBurstRaster: BurstIndex %d out of range [1, %d]', bi, numel(peak_locs));
    peak_locs = peak_locs(bi);
end

n_bursts = numel(peak_locs);
n_units  = numel(ud_indices);

if n_bursts == 0
    warning('EIAnalyzer:PlotBurstRaster', ...
        'No bursts remain after filtering — skipping culture %d.', culture_idx);
    ax = gobjects(0);
    return
end

% Convert peak bin indices to seconds using Activity time axis
peak_times = act.x(peak_locs);  % (1 x n_bursts) seconds

% ── Sort units ────────────────────────────────────────────────────────────
% Sort key: (1) cell type group, (2) mean spike count per burst (descending)
mean_count = zeros(n_units, 1);
for u = 1:n_units
    sp = ud(ud_indices(u)).SpikeTimes;
    cnt = 0;
    for b = 1:n_bursts
        t0  = peak_times(b) - half_w;
        t1  = peak_times(b) + half_w;
        cnt = cnt + sum(sp >= t0 & sp <= t1);
    end
    mean_count(u) = cnt / n_bursts;
end

type_group = zeros(n_units, 1);
type_group(labels == 1)   = 1;   % excitatory first
type_group(labels == 2)   = 2;   % inhibitory second
type_group(isnan(labels)) = 3;   % unclassified last
[~, sort_idx] = sortrows([type_group, -mean_count]);

sorted_ud_idx = ud_indices(sort_idx);
sorted_labels = labels(sort_idx);
n_exc  = sum(labels == 1);
n_inh  = sum(labels == 2);

% ── Color per unit ────────────────────────────────────────────────────────
exc_color = [0.12, 0.47, 0.71];  % blue
inh_color = [0.84, 0.19, 0.15];  % red
unc_color = [0.60, 0.60, 0.60];  % grey

unit_colors = repmat(unc_color, n_units, 1);
unit_colors(sorted_labels == 1, :) = repmat(exc_color, sum(sorted_labels == 1), 1);
unit_colors(sorted_labels == 2, :) = repmat(inh_color, sum(sorted_labels == 2), 1);

% ── Axes ─────────────────────────────────────────────────────────────────
if isempty(opts.Parent)
    figure('Color', 'w');
    ax = axes();
else
    ax = opts.Parent;
end

hold(ax, 'on');

% Alpha: spread across bursts so individual bursts are visible; min 0.15
pt_alpha = max(0.15, min(1, 2 / max(n_bursts, 1)));

% ── Plot raster ───────────────────────────────────────────────────────────
for u = 1:n_units
    sp  = ud(sorted_ud_idx(u)).SpikeTimes;
    col = unit_colors(u, :);
    for b = 1:n_bursts
        t0   = peak_times(b) - half_w;
        t1   = peak_times(b) + half_w;
        sp_b = sp(sp >= t0 & sp <= t1) - peak_times(b);  % relative time (s)
        if ~isempty(sp_b)
            scatter(ax, sp_b, repmat(u, size(sp_b)), 2, col, 'filled', ...
                'MarkerFaceAlpha', pt_alpha, 'MarkerEdgeColor', 'none');
        end
    end
end

% Mark burst peak with dashed vertical line
xline(ax, 0, '--k', 'Alpha', 0.35, 'LineWidth', 0.8);

% Horizontal separator between E and I populations
if n_exc > 0 && n_inh > 0
    yline(ax, n_exc + 0.5, '-', 'Color', [0.4 0.4 0.4], 'Alpha', 0.5, 'LineWidth', 0.8);
end

xlim(ax, [-half_w, half_w]);
ylim(ax, [0.5, n_units + 0.5]);
set(ax, 'YDir', 'reverse');  % unit 1 at top
box(ax, 'off');

xlabel(ax, 'Time relative to burst peak (s)');
ylabel(ax, 'Unit (sorted: E → I → unclassified)');

if isequal(opts.BurstIndex, 'all')
    burst_desc = sprintf('all %d accepted bursts', n_bursts);
else
    burst_desc = sprintf('burst %d', opts.BurstIndex);
end
title(ax, sprintf('Culture %d — burst raster (%s, %d E, %d I)', ...
    culture_idx, burst_desc, n_exc, n_inh), 'FontWeight', 'normal');

setAllFontSizes(ancestor(ax, 'figure'), eia.Parameters.Plotting.FontSize);
end

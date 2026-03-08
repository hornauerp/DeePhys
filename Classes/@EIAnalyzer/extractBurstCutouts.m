function extractBurstCutouts(eia)
% EXTRACTBURSTCUTOUTS  Extract and align burst cutouts around peak network activity.
%
% For each culture:
%   1. Detects burst peaks in excitatory-dominant activity
%   2. Extracts ±PeakCutout bins around each peak
%   3. Aligns all bursts to the mean burst shape via cross-correlation
%   4. Clusters bursts into 2 types using t-SNE + kmedoids
%
% Returns four variants (stored in a cell array indexed by SelectedVariant):
%   1 — all bursts (cross-correlation aligned)
%   2 — t-SNE cluster 1 only
%   3 — t-SNE cluster 2 only
%   4 — both clusters aligned to each other (cluster 2 shifted to match cluster 1)
%
% Requires eia.Activity and eia.BurstState (run computeActivity + detectBursts first).
% Sets: eia.BurstCutouts (selected variant per culture)

arguments
    eia EIAnalyzer
end

assert(~isempty(eia.BurstState), 'Run detectBursts() before extractBurstCutouts()');

p       = eia.Parameters.BurstDetection;
sel     = p.SelectedVariant;
cutouts = repmat(struct('total',[],'i_mat',[],'e_mat',[]), 1, length(eia.Activity));

for c = 1:length(eia.Activity)
    act = eia.Activity(c);
    z   = eia.BurstState{c};

    % Detect peaks in burst-state excitatory activity
    only_burst = act.total .* (z' == 2);
    [~, locs]  = findpeaks(only_burst, ...
        'MinPeakHeight',     p.PeakHeight, ...
        'MinPeakProminence', p.PeakProminence, ...
        'MinPeakDistance',   6);

    % Remove peaks too close to the edges
    locs(locs <= p.PeakCutout | locs > length(only_burst) - p.PeakCutout) = [];

    n_bursts = length(locs);
    if n_bursts < 4
        % Not enough bursts to cluster; leave empty
        continue
    end

    w         = 2 * p.PeakCutout + 1;
    burst_mat = nan(n_bursts, w);
    i_mat     = nan(n_bursts, w);
    e_mat     = nan(n_bursts, w);

    for l = 1:n_bursts
        idx           = locs(l) - p.PeakCutout : locs(l) + p.PeakCutout;
        burst_mat(l,:) = act.total(idx);
        i_mat(l,:)    = act.I(idx);
        e_mat(l,:)    = act.E(idx);
    end

    % Drop bursts with NaN bins
    valid     = ~any(isnan(burst_mat), 2);
    burst_mat = burst_mat(valid,:);
    i_mat     = i_mat(valid,:);
    e_mat     = e_mat(valid,:);

    if size(burst_mat, 1) < 4; continue; end

    % Cross-correlation alignment to mean burst shape
    norm_burst = burst_mat ./ max(burst_mat, [], 2);
    mean_burst = mean(norm_burst, 1);
    for b = 1:size(burst_mat, 1)
        [~, max_lag]  = max(xcorr(norm_burst(b,:), mean_burst, p.PeakCutout));
        shift         = (p.PeakCutout + 1) - max_lag;
        burst_mat(b,:) = circshift(burst_mat(b,:), shift);
        i_mat(b,:)    = circshift(i_mat(b,:),     shift);
        e_mat(b,:)    = circshift(e_mat(b,:),     shift);
    end

    % Build four burst variants
    ba = cell(1, 4);
    ba{1} = struct('total', burst_mat, 'i_mat', i_mat, 'e_mat', e_mat);

    try
        ts  = tsne(burst_mat, 'Standardize', true, 'NumDimensions', 5);
        idx = kmedoids(ts, 2);

        ba{2} = struct('total', burst_mat(idx==1,:), 'i_mat', i_mat(idx==1,:), 'e_mat', e_mat(idx==1,:));
        ba{3} = struct('total', burst_mat(idx==2,:), 'i_mat', i_mat(idx==2,:), 'e_mat', e_mat(idx==2,:));

        % Align cluster 2 to cluster 1 by cross-correlation
        [~, align_lag] = max(xcorr(mean(ba{3}.total, 1), mean(ba{2}.total, 1), p.PeakCutout));
        shift2         = (p.PeakCutout + 1) - align_lag;
        bm4 = burst_mat; im4 = i_mat; em4 = e_mat;
        id2 = find(idx == 2);
        for b = 1:length(id2)
            bm4(id2(b),:) = circshift(burst_mat(id2(b),:), shift2);
            im4(id2(b),:) = circshift(i_mat(id2(b),:),     shift2);
            em4(id2(b),:) = circshift(e_mat(id2(b),:),     shift2);
        end
        ba{4} = struct('total', bm4, 'i_mat', im4, 'e_mat', em4);
    catch
        % Clustering failed — use all-bursts variant for all slots
        [ba{2}, ba{3}, ba{4}] = deal(ba{1});
    end

    cutouts(c) = ba{sel};
end

eia.BurstCutouts = cutouts;
end

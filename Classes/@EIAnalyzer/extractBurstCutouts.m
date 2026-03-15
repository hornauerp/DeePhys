function extractBurstCutouts(eia)
% EXTRACTBURSTCUTOUTS  Extract and align burst cutouts around peak network activity.
%
% For each culture:
%   1. Detects burst peaks in total activity (optionally masked to burst-state bins)
%   2. Extracts +/-PeakCutout bins around each peak
%   3. Aligns all bursts to the mean burst shape via cross-correlation
%   4. Filters bursts using the selected method:
%        "correlation" — iterative correlation-based filtering (deterministic)
%        "tsne"        — t-SNE + kmedoids clustering (seeded for reproducibility)
%        "none"        — keep all aligned bursts
%
% Variant structure (indexed by SelectedVariant):
%   For "tsne":
%     1 — all bursts (xcorr aligned)
%     2 — t-SNE cluster 1 only
%     3 — t-SNE cluster 2 only
%     4 — both clusters aligned to each other
%   For "correlation":
%     1 — all bursts (xcorr aligned, before filtering)
%     2 — rejected bursts (below correlation threshold)
%     3 — accepted bursts (above correlation threshold)
%     4 — accepted bursts (same as 3, canonical output)
%   For "none":
%     1–4 — all bursts (xcorr aligned)
%
% Parameters (from eia.Parameters.BurstDetection):
%   MaskPeaks            — true: peak detection only in burst-state bins (z==2)
%                          false: peak detection on full activity trace (main-branch default)
%   BurstFiltering       — "correlation" | "tsne" | "none"
%   CorrelationThreshold — min Pearson r with mean template (for "correlation")
%   MaxFilterIter        — max refinement iterations (for "correlation")
%   RNGSeed              — seed for t-SNE + kmedoids reproducibility
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

    % ── Peak detection ───────────────────────────────────────────────────
    if p.MaskPeaks
        peak_trace = act.total .* (z == 2);
    else
        peak_trace = act.total;
    end

    [~, locs] = findpeaks(peak_trace, ...
        'MinPeakHeight',     p.PeakHeight, ...
        'MinPeakProminence', p.PeakProminence, ...
        'MinPeakDistance',   p.MinPeakDistance);

    % Remove peaks too close to the edges
    locs(locs <= p.PeakCutout | locs > length(peak_trace) - p.PeakCutout) = [];

    n_bursts = length(locs);
    if n_bursts < 4; continue; end

    % ── Extract cutouts ──────────────────────────────────────────────────
    w         = 2 * p.PeakCutout + 1;
    burst_mat = nan(n_bursts, w);
    i_mat     = nan(n_bursts, w);
    e_mat     = nan(n_bursts, w);

    for l = 1:n_bursts
        idx            = locs(l) - p.PeakCutout : locs(l) + p.PeakCutout;
        burst_mat(l,:) = act.total(idx);
        i_mat(l,:)     = act.I(idx);
        e_mat(l,:)     = act.E(idx);
    end

    % Drop bursts with NaN bins
    valid     = ~any(isnan(burst_mat), 2);
    burst_mat = burst_mat(valid,:);
    i_mat     = i_mat(valid,:);
    e_mat     = e_mat(valid,:);

    if size(burst_mat, 1) < 4; continue; end

    % ── Cross-correlation alignment to mean burst shape ──────────────────
    [burst_mat, i_mat, e_mat] = xcorrAlign(burst_mat, i_mat, e_mat, p.PeakCutout);

    % ── Build variants by filtering method ───────────────────────────────
    ba = cell(1, 4);
    ba{1} = struct('total', burst_mat, 'i_mat', i_mat, 'e_mat', e_mat);

    switch p.BurstFiltering
        case "tsne"
            ba = filterByTSNE(ba, burst_mat, i_mat, e_mat, p);

        case "correlation"
            ba = filterByCorrelation(ba, burst_mat, i_mat, e_mat, p);

        case "none"
            [ba{2}, ba{3}, ba{4}] = deal(ba{1});

        otherwise
            warning('EIAnalyzer:extractBurstCutouts', ...
                'Unknown BurstFiltering method "%s" — keeping all bursts.', p.BurstFiltering);
            [ba{2}, ba{3}, ba{4}] = deal(ba{1});
    end

    cutouts(c) = ba{sel};
end

eia.BurstCutouts = cutouts;
end

%% ── Helper functions ────────────────────────────────────────────────────────

function [burst_mat, i_mat, e_mat] = xcorrAlign(burst_mat, i_mat, e_mat, peak_cutout)
% Align bursts to population mean via cross-correlation (NaN-padded shifts).
    norm_burst = burst_mat ./ max(burst_mat, [], 2);
    mean_burst = mean(norm_burst, 1, 'omitnan');
    for b = 1:size(burst_mat, 1)
        [~, max_lag]   = max(xcorr(norm_burst(b,:), mean_burst, peak_cutout));
        shift          = (peak_cutout + 1) - max_lag;
        burst_mat(b,:) = nanShift(burst_mat(b,:), shift);
        i_mat(b,:)     = nanShift(i_mat(b,:),     shift);
        e_mat(b,:)     = nanShift(e_mat(b,:),     shift);
    end
end

function ba = filterByTSNE(ba, burst_mat, i_mat, e_mat, p)
% t-SNE + kmedoids clustering into 2 burst types (seeded for reproducibility).
    try
        rng(p.RNGSeed);
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
            bm4(id2(b),:) = nanShift(burst_mat(id2(b),:), shift2);
            im4(id2(b),:) = nanShift(i_mat(id2(b),:),     shift2);
            em4(id2(b),:) = nanShift(e_mat(id2(b),:),     shift2);
        end
        ba{4} = struct('total', bm4, 'i_mat', im4, 'e_mat', em4);
    catch ME
        warning('EIAnalyzer:filterByTSNE', 't-SNE/kmedoids filtering failed: %s. Keeping all bursts.', ME.message);
        [ba{2}, ba{3}, ba{4}] = deal(ba{1});
    end
end

function ba = filterByCorrelation(ba, burst_mat, i_mat, e_mat, p)
% Iterative correlation-based filtering: keep bursts that correlate well
% with the population mean template.  Deterministic and robust.
%
%   1. Compute mean burst template (normalised)
%   2. Correlate each burst with template
%   3. Remove bursts below CorrelationThreshold
%   4. Recompute template from surviving bursts
%   5. Repeat until stable or MaxFilterIter reached
%
% Variants:
%   2 — rejected bursts
%   3 — accepted bursts
%   4 — accepted bursts (canonical output)

    keep = true(size(burst_mat, 1), 1);

    for iter = 1:p.MaxFilterIter
        template = mean(burst_mat(keep,:), 1, 'omitnan');
        norm_template = template ./ max(abs(template));
        if all(norm_template == 0); break; end

        r = nan(size(burst_mat, 1), 1);
        for b = 1:size(burst_mat, 1)
            row = burst_mat(b,:);
            valid = ~isnan(row) & ~isnan(template);
            if sum(valid) < 10; r(b) = 0; continue; end
            cc = corrcoef(row(valid), template(valid));
            r(b) = cc(1,2);
        end

        new_keep = r >= p.CorrelationThreshold;
        if isequal(new_keep, keep); break; end  % converged
        if sum(new_keep) < 4
            % Too few bursts survive — relax to keep at least 4
            [~, rank] = sort(r, 'descend');
            new_keep = false(size(keep));
            new_keep(rank(1:min(4, end))) = true;
            keep = new_keep;
            break
        end
        keep = new_keep;
    end

    reject = ~keep;
    ba{2} = struct('total', burst_mat(reject,:), 'i_mat', i_mat(reject,:), 'e_mat', e_mat(reject,:));
    ba{3} = struct('total', burst_mat(keep,:),   'i_mat', i_mat(keep,:),   'e_mat', e_mat(keep,:));
    ba{4} = struct('total', burst_mat(keep,:),   'i_mat', i_mat(keep,:),   'e_mat', e_mat(keep,:));
end

function out = nanShift(row, shift)
% NANSHIFT  Shift a row vector by shift bins, NaN-padding the vacated edge.
    n   = length(row);
    out = nan(1, n);
    if shift == 0
        out = row;
    elseif shift > 0
        out(shift+1:end) = row(1:end-shift);
    else  % shift < 0
        out(1:end+shift) = row(1-shift:end);
    end
end

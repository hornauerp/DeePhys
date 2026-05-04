function [network_table, unit_table] = computeSpatialFeatures( ...
    electrode_coordinates, units, cell_type_labels, unit_feature_table, params)
% COMPUTESPATIALFEATURES  Spatial spread, E/I mixing, and firing rate gradient features.
%
%   [network_table, unit_table] = computeSpatialFeatures( ...
%       electrode_coordinates, units, cell_type_labels, unit_feature_table, params)
%
% Computes three groups of features:
%   A. Spatial spread & coverage  — convex hull, pairwise distances, centroid spread
%   B. Nearest-neighbor E/I composition — local mixing of cell types
%   E. Firing rate gradient       — spatial autocorrelation, center-periphery ratio
%
% INPUTS:
%   electrode_coordinates - (N_channels x 2) XY electrode positions in micrometers
%   units                 - (1 x N_units) UnitData array with .ReferenceElectrode
%   cell_type_labels      - (1 x N_units) 1=excitatory, 2=inhibitory, NaN=unclassified
%                           (may be empty or all-NaN — cell-type features skipped)
%   unit_feature_table    - table with FiringRate column (optional; pass [] to skip FR features)
%   params                - struct with optional fields:
%                             .NNk            (default 5)  k for nearest-neighbour composition
%                             .RipleyEnable   (default false) run Ripley's K/L
%                             .RipleyMaxRadius (default 500 um)
%                             .RipleyNBins    (default 50)
%                             .RipleyNSurrogates (default 100)
%
% OUTPUTS:
%   network_table - (1 x F) table of network-level spatial features
%   unit_table    - (N_units x F) table of per-unit spatial features

    arguments
        electrode_coordinates (:,2) double
        units                 (1,:)
        cell_type_labels      (1,:) double     = []
        unit_feature_table              = []
        params                struct   = struct()
    end

    % --- Defaults ---
    p = parseStructParameters(struct( ...
        'NNk',                    5, ...
        'RipleyEnable',           false, ...
        'RipleyMaxRadius',        500, ...
        'RipleyNBins',            50, ...
        'RipleyNSurrogates',      100, ...
        'KernelBandwidth',        [], ...
        'BurstMinUnitsForWavefront', []), params);

    n_units = numel(units);
    network_table = table();
    unit_table    = table();

    if n_units < 2
        unit_table = makeNaNUnitTable(n_units);
        return
    end

    % Gather unit XY positions
    ref_els = arrayfun(@(u) u.ReferenceElectrode, units);
    xy      = electrode_coordinates(ref_els, :);   % (N_units x 2)

    has_labels = ~isempty(cell_type_labels) && numel(cell_type_labels) == n_units ...
                 && any(~isnan(cell_type_labels));
    if has_labels
        exc_idx = find(cell_type_labels == 1);
        inh_idx = find(cell_type_labels == 2);
    else
        exc_idx = [];
        inh_idx = [];
    end

    % =========================================================
    % A. Spatial spread & coverage
    % =========================================================
    nw_spread = computeSpread(xy, exc_idx, inh_idx, electrode_coordinates);
    network_table = [network_table, struct2table(nw_spread)];

    % =========================================================
    % B. Nearest-neighbour cell type composition
    % =========================================================
    if has_labels && n_units >= p.NNk + 1
        [nw_nn, unit_nn] = computeNNComposition(xy, cell_type_labels, exc_idx, inh_idx, p.NNk);
        network_table = [network_table, struct2table(nw_nn)];
        unit_table    = [unit_table, unit_nn];
    else
        nw_nn = makeNaNNNStruct();
        network_table = [network_table, struct2table(nw_nn)];
        unit_table    = [unit_table, table( ...
            NaN(n_units,1), NaN(n_units,1), ...
            'VariableNames', {'FractionExcNN', 'FractionInhNN'})];
    end

    % =========================================================
    % E. Firing rate gradient & centroid distance
    % =========================================================
    centroid      = mean(xy, 1);
    dist_centroid = sqrt(sum((xy - centroid).^2, 2));
    unit_table    = [unit_table, table(dist_centroid, 'VariableNames', {'DistFromCentroid'})];

    if ~isempty(unit_feature_table) && ...
            ismember('FiringRate', unit_feature_table.Properties.VariableNames)
        fr = unit_feature_table.FiringRate;
        nw_fr = computeFRGradient(fr, xy, dist_centroid, exc_idx, inh_idx);
        network_table = [network_table, struct2table(nw_fr)];
    end

    % =========================================================
    % G. Ripley's K/L (optional, expensive)
    % =========================================================
    if p.RipleyEnable && n_units >= 10
        study_area = convhullArea(xy);
        [L_max_all, scale_all] = computeRipleysK(xy, study_area, p);
        rip_nw.RipleysL_max    = L_max_all;
        rip_nw.ClusterScale    = scale_all;

        if numel(exc_idx) >= 10
            [L_max_e, scale_e] = computeRipleysK(xy(exc_idx,:), study_area, p);
            rip_nw.RipleysL_max_E = L_max_e;
            rip_nw.ClusterScale_E = scale_e;
        else
            rip_nw.RipleysL_max_E = NaN;
            rip_nw.ClusterScale_E = NaN;
        end

        if numel(inh_idx) >= 10
            [L_max_i, scale_i] = computeRipleysK(xy(inh_idx,:), study_area, p);
            rip_nw.RipleysL_max_I = L_max_i;
            rip_nw.ClusterScale_I = scale_i;
        else
            rip_nw.RipleysL_max_I = NaN;
            rip_nw.ClusterScale_I = NaN;
        end

        network_table = [network_table, struct2table(rip_nw)];
    end
end


% =========================================================
% Spatial spread helpers
% =========================================================
function nw = computeSpread(xy, exc_idx, inh_idx, all_electrodes)
    nw.ConvexHullArea = convhullArea(xy);
    nw.ChipCoverage   = nw.ConvexHullArea / max(convhullArea(all_electrodes), eps);

    D = squareform(pdist(xy));
    ut = D(triu(true(size(D)), 1));
    nw.MeanPairwiseDistance = mean(ut, 'omitnan');

    centroid = mean(xy, 1);
    nw.CentroidSpread = mean(sqrt(sum((xy - centroid).^2, 2)), 'omitnan');

    % Type-stratified
    nw.ConvexHullArea_E         = typeHullArea(xy, exc_idx);
    nw.ConvexHullArea_I         = typeHullArea(xy, inh_idx);
    nw.MeanPairwiseDistance_EE  = typeMeanDist(D, exc_idx, exc_idx, true);
    nw.MeanPairwiseDistance_II  = typeMeanDist(D, inh_idx, inh_idx, true);
    nw.MeanPairwiseDistance_EI  = typeMeanDist(D, exc_idx, inh_idx, false);
    nw.CentroidSpread_E         = typeCentroidSpread(xy, exc_idx);
    nw.CentroidSpread_I         = typeCentroidSpread(xy, inh_idx);
end

function a = convhullArea(xy)
    if size(xy, 1) < 3
        a = NaN;
        return
    end
    try
        k = convhull(xy(:,1), xy(:,2));
        a = polyarea(xy(k,1), xy(k,2));
    catch
        a = NaN;
    end
end

function a = typeHullArea(xy, idx)
    if numel(idx) < 3
        a = NaN;
    else
        a = convhullArea(xy(idx,:));
    end
end

function m = typeMeanDist(D, row_idx, col_idx, upper_tri_only)
    if isempty(row_idx) || isempty(col_idx)
        m = NaN;
        return
    end
    sub = D(row_idx, col_idx);
    if upper_tri_only && isequal(row_idx(:), col_idx(:))
        sub = sub(triu(true(size(sub)), 1));
    else
        sub = sub(:);
    end
    m = mean(sub, 'omitnan');
end

function s = typeCentroidSpread(xy, idx)
    if numel(idx) < 2
        s = NaN;
        return
    end
    sub = xy(idx,:);
    c   = mean(sub, 1);
    s   = mean(sqrt(sum((sub - c).^2, 2)), 'omitnan');
end


% =========================================================
% Nearest-neighbour composition
% =========================================================
function [nw, unit_tbl] = computeNNComposition(xy, labels, exc_idx, inh_idx, k)
    n = size(xy, 1);
    frac_exc = NaN(n, 1);
    frac_inh = NaN(n, 1);

    % Build distance matrix once
    D = squareform(pdist(xy));
    D(1:n+1:end) = Inf;  % exclude self

    for u = 1:n
        [~, nn_order] = sort(D(u,:));
        nn_idx = nn_order(1:k);
        nn_labels = labels(nn_idx);
        valid = ~isnan(nn_labels);
        if sum(valid) == 0, continue; end
        frac_exc(u) = sum(nn_labels(valid) == 1) / sum(valid);
        frac_inh(u) = sum(nn_labels(valid) == 2) / sum(valid);
    end

    unit_tbl = table(frac_exc, frac_inh, ...
        'VariableNames', {'FractionExcNN', 'FractionInhNN'});

    nw.SpatialMixingIndex   = NaN;
    nw.MeanFractionExcNN_E  = NaN;
    nw.MeanFractionExcNN_I  = NaN;

    labeled = ~isnan(labels);
    if sum(labeled) > 0
        % Mixing index: mean fraction of unlike-type neighbors (for labeled units)
        like_frac = NaN(n, 1);
        for u = 1:n
            if isnan(labels(u)), continue; end
            if labels(u) == 1
                like_frac(u) = frac_exc(u);
            else
                like_frac(u) = frac_inh(u);
            end
        end
        unlike_frac = 1 - like_frac(labeled);
        nw.SpatialMixingIndex = mean(unlike_frac, 'omitnan');
    end

    if ~isempty(exc_idx)
        nw.MeanFractionExcNN_E = mean(frac_exc(exc_idx), 'omitnan');
    end
    if ~isempty(inh_idx)
        nw.MeanFractionExcNN_I = mean(frac_exc(inh_idx), 'omitnan');
    end
end

function nw = makeNaNNNStruct()
    nw.SpatialMixingIndex  = NaN;
    nw.MeanFractionExcNN_E = NaN;
    nw.MeanFractionExcNN_I = NaN;
end


% =========================================================
% Firing rate gradient
% =========================================================
function nw = computeFRGradient(fr, xy, dist_centroid, exc_idx, inh_idx)
    n = numel(fr);

    % Moran's I using inverse-distance weights
    D = squareform(pdist(xy));
    D(1:n+1:end) = Inf;
    W = 1 ./ D;
    W(isinf(W)) = 0;
    nw.SpatialFRMoransI = moransI(fr, W);

    % Center-periphery split (median distance)
    med_dist = median(dist_centroid);
    center_mask = dist_centroid <= med_dist;
    periph_mask = dist_centroid >  med_dist;

    center_fr = mean(fr(center_mask), 'omitnan');
    periph_fr = mean(fr(periph_mask), 'omitnan');
    nw.CenterPeripheryFR_ratio = safeFRRatio(center_fr, periph_fr);

    nw.CenterPeripheryFR_ratio_E = NaN;
    nw.CenterPeripheryFR_ratio_I = NaN;
    if ~isempty(exc_idx)
        c_e = mean(fr(intersect(exc_idx, find(center_mask))), 'omitnan');
        p_e = mean(fr(intersect(exc_idx, find(periph_mask))), 'omitnan');
        nw.CenterPeripheryFR_ratio_E = safeFRRatio(c_e, p_e);
    end
    if ~isempty(inh_idx)
        c_i = mean(fr(intersect(inh_idx, find(center_mask))), 'omitnan');
        p_i = mean(fr(intersect(inh_idx, find(periph_mask))), 'omitnan');
        nw.CenterPeripheryFR_ratio_I = safeFRRatio(c_i, p_i);
    end

    % DistFromCentroid added by caller
end

function I = moransI(x, W)
% Moran's I spatial autocorrelation.
    x = x(:);
    valid = ~isnan(x);
    if sum(valid) < 4
        I = NaN;
        return
    end
    Wv = W(valid, valid);
    xv = x(valid);
    n  = sum(valid);
    xm = mean(xv);
    z  = xv - xm;
    S0 = sum(Wv(:));
    if S0 == 0
        I = NaN;
        return
    end
    I  = (n / S0) * (z' * (Wv * z)) / (z' * z);
end


% =========================================================
% Utilities
% =========================================================
function tbl = makeNaNUnitTable(n)
    tbl = table( ...
        NaN(n,1), NaN(n,1), NaN(n,1), ...
        'VariableNames', {'FractionExcNN', 'FractionInhNN', 'DistFromCentroid'});
end

function r = safeFRRatio(numer, denom)
% Return NaN when denominator is zero or NaN to avoid eps-inflated ratios.
    if isnan(denom) || denom == 0
        r = NaN;
    else
        r = numer / denom;
    end
end

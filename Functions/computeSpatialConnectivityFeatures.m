function [network_table, unit_table] = computeSpatialConnectivityFeatures( ...
    electrode_coordinates, units, connectivity, cell_type_labels, params)
% COMPUTESPATIALCONNECTIVITYFEATURES  Distance-dependent functional connectivity.
%
%   [network_table, unit_table] = computeSpatialConnectivityFeatures( ...
%       electrode_coordinates, units, connectivity, cell_type_labels, params)
%
% Correlates STTC weighted connectivity with Euclidean distance between units.
% Fits an exponential decay model and computes distance-FC correlation per
% cell-type pair category (EE, II, EI).
%
% INPUTS:
%   electrode_coordinates - (N_channels x 2) XY positions in micrometers
%   units                 - (1 x N_units) UnitData array with .ReferenceElectrode
%   connectivity          - struct with .STTC.wu (N x N weighted undirected matrix)
%   cell_type_labels      - (1 x N_units) 1=exc, 2=inh, NaN=unclassified (may be [])
%   params                - struct (optional); field .NNk for per-unit near-neighbor FC
%
% OUTPUTS:
%   network_table - (1 x F) table: SpatialDecayTau*, DistanceFCCorrelation*
%   unit_table    - (N x 1) table: MeanFC_NearNeighbors

    arguments
        electrode_coordinates (:,2) double
        units                 (1,:)
        connectivity          struct
        cell_type_labels      (1,:) double = []
        params                struct = struct()
    end

    p = parseStructParameters(struct( ...
        'NNk',                       5, ...
        'RipleyEnable',              [], ...
        'RipleyMaxRadius',           [], ...
        'RipleyNBins',               [], ...
        'RipleyNSurrogates',         [], ...
        'KernelBandwidth',           [], ...
        'BurstMinUnitsForWavefront', []), params);

    n_units = numel(units);
    nw = makeNaNStruct();
    unit_table = table(NaN(n_units,1), 'VariableNames', {'MeanFC_NearNeighbors'});
    network_table = struct2table(nw);

    if ~isfield(connectivity, 'STTC') || ~isfield(connectivity.STTC, 'wu')
        return
    end
    W = connectivity.STTC.wu;
    if isempty(W) || size(W,1) ~= n_units
        return
    end

    ref_els = arrayfun(@(u) u.ReferenceElectrode, units);
    xy      = electrode_coordinates(ref_els, :);

    D = squareform(pdist(xy));   % (N x N) pairwise distances in um

    % Upper-triangle pairs only (symmetric matrix)
    ut_mask = triu(true(n_units), 1);
    d_vec   = D(ut_mask);
    w_vec   = W(ut_mask);
    valid   = ~isnan(w_vec) & ~isnan(d_vec) & d_vec > 0;

    if sum(valid) >= 5
        nw.DistanceFCCorrelation = corr(d_vec(valid), w_vec(valid), 'type', 'Pearson');
        nw.SpatialDecayTau       = fitExpDecay(d_vec(valid), w_vec(valid));
    end

    % Cell-type stratified decay
    if ~isempty(cell_type_labels) && numel(cell_type_labels) == n_units
        exc_idx = find(cell_type_labels == 1);
        inh_idx = find(cell_type_labels == 2);

        nw.SpatialDecayTau_EE = pairDecay(D, W, exc_idx, exc_idx, true);
        nw.SpatialDecayTau_II = pairDecay(D, W, inh_idx, inh_idx, true);
        nw.SpatialDecayTau_EI = pairDecay(D, W, exc_idx, inh_idx, false);
    end

    % Per-unit: mean STTC with k nearest spatial neighbors
    D_self_inf = D;
    D_self_inf(1:n_units+1:end) = Inf;
    mean_fc_nn = NaN(n_units, 1);
    k = min(p.NNk, n_units - 1);
    for u = 1:n_units
        [~, nn_order] = sort(D_self_inf(u,:));
        nn_idx = nn_order(1:k);
        mean_fc_nn(u) = mean(W(u, nn_idx), 'omitnan');
    end
    unit_table = table(mean_fc_nn, 'VariableNames', {'MeanFC_NearNeighbors'});

    network_table = struct2table(nw);
end


function tau = pairDecay(D, W, row_idx, col_idx, upper_tri_only)
    if isempty(row_idx) || isempty(col_idx)
        tau = NaN;
        return
    end
    d_sub = D(row_idx, col_idx);
    w_sub = W(row_idx, col_idx);
    if upper_tri_only && isequal(row_idx(:), col_idx(:))
        mask  = triu(true(numel(row_idx)), 1);
        d_vec = d_sub(mask);
        w_vec = w_sub(mask);
    else
        d_vec = d_sub(:);
        w_vec = w_sub(:);
    end
    valid = ~isnan(w_vec) & ~isnan(d_vec) & d_vec > 0;
    if sum(valid) < 5
        tau = NaN;
        return
    end
    tau = fitExpDecay(d_vec(valid), w_vec(valid));
end

function tau = fitExpDecay(d, w)
% Least-squares fit of w = a*exp(-d/tau) + c over all pairs (zero and negative w included:
% a fit on the positive values only is biased towards slow decay). For a given tau, a and
% c are linear, so tau is found by a profile search between the smallest and largest
% distance. NaN when the coupling does not decay (a <= 0), the decay explains no more than
% a constant (F-test, p >= 0.05), or tau reaches either end of the measured distances
% (decay not resolved within the measured range).
    d = d(:); w = w(:);
    tau = NaN;
    if numel(d) < 4 || range(d) == 0, return; end
    lo = log(max(min(d), eps)); hi = log(max(d));
    sse = @(lt) profileSSE(lt, d, w);
    grid = linspace(lo, hi, 60);
    s = arrayfun(sse, grid);
    [~, k] = min(s);
    lt = fminbnd(sse, grid(max(k - 1, 1)), grid(min(k + 1, end)));
    tol = 1e-3 * (hi - lo);
    if lt >= hi - tol || lt <= lo + tol, return; end
    [s1, a] = profileSSE(lt, d, w);
    % The decay must explain more than a constant (nested F-test, 2 extra parameters)
    n = numel(w); s0 = sum((w - mean(w)).^2);
    if n <= 3 || s1 <= 0, return; end
    p = 1 - fcdf(((s0 - s1) / 2) / (s1 / (n - 3)), 2, n - 3);
    if a > 0 && p < 0.05, tau = exp(lt); end
end

function [s, a] = profileSSE(log_tau, d, w)
    X = [exp(-d / exp(log_tau)), ones(size(d))];
    b = X \ w;
    s = sum((w - X * b).^2);
    a = b(1);
end

function nw = makeNaNStruct()
    nw.SpatialDecayTau       = NaN;
    nw.SpatialDecayTau_EE    = NaN;
    nw.SpatialDecayTau_II    = NaN;
    nw.SpatialDecayTau_EI    = NaN;
    nw.DistanceFCCorrelation = NaN;
end

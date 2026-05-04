function [network_table, viz_data] = computeSpatialEIBalance( ...
    electrode_coordinates, units, cell_type_labels, params)
% COMPUTESPATIALEIBALANCE  Local E/I balance map and spatial heterogeneity metrics.
%
%   [network_table, viz_data] = computeSpatialEIBalance( ...
%       electrode_coordinates, units, cell_type_labels, params)
%
% Computes a kernel-smoothed estimate of the local excitatory fraction at
% each grid point on the chip. Derives summary statistics about spatial
% variability of E/I balance and its spatial autocorrelation.
%
% INPUTS:
%   electrode_coordinates - (N_channels x 2) XY positions in micrometers
%   units                 - (1 x N_units) UnitData array with .ReferenceElectrode
%   cell_type_labels      - (1 x N_units) 1=exc, 2=inh, NaN=unclassified
%   params                - struct with optional fields:
%                             .KernelBandwidth (default 100 um)
%                             .GridSpacing     (default 50 um)
%
% OUTPUTS:
%   network_table - (1 x 2) table: SpatialEIVariability, SpatialEIMoransI
%   viz_data      - struct with fields:
%                     .grid_x, .grid_y  — meshgrid coordinates
%                     .exc_fraction     — local excitatory fraction at each grid point
%                     .unit_xy          — unit positions
%                     .unit_labels      — cell type labels

    arguments
        electrode_coordinates (:,2) double
        units                 (1,:)
        cell_type_labels      (1,:) double
        params                struct = struct()
    end

    p = parseStructParameters(struct( ...
        'NNk',                       [], ...
        'RipleyEnable',              [], ...
        'RipleyMaxRadius',           [], ...
        'RipleyNBins',               [], ...
        'RipleyNSurrogates',         [], ...
        'KernelBandwidth',           100, ...
        'GridSpacing',               50, ...
        'BurstMinUnitsForWavefront', []), params);

    nw.SpatialEIVariability = NaN;
    nw.SpatialEIMoransI     = NaN;
    network_table = struct2table(nw);
    viz_data = struct();

    ref_els = arrayfun(@(u) u.ReferenceElectrode, units);
    xy      = electrode_coordinates(ref_els, :);
    labels  = cell_type_labels(:);

    exc_mask = labels == 1;
    inh_mask = labels == 2;
    if sum(exc_mask) < 2 || sum(inh_mask) < 2
        return
    end

    % Build evaluation grid
    x_range = [min(xy(:,1)), max(xy(:,1))];
    y_range = [min(xy(:,2)), max(xy(:,2))];
    gx = x_range(1):p.GridSpacing:x_range(2);
    gy = y_range(1):p.GridSpacing:y_range(2);
    if numel(gx) < 2 || numel(gy) < 2
        return
    end
    [GX, GY] = meshgrid(gx, gy);
    grid_pts  = [GX(:), GY(:)];
    n_grid    = size(grid_pts, 1);

    h = p.KernelBandwidth;
    exc_density = zeros(n_grid, 1);
    inh_density = zeros(n_grid, 1);

    for g = 1:n_grid
        gp = grid_pts(g,:);
        d2 = sum((xy - gp).^2, 2);
        w  = exp(-d2 / (2 * h^2));
        exc_density(g) = sum(w .* exc_mask);
        inh_density(g) = sum(w .* inh_mask);
    end

    total_density = exc_density + inh_density;
    valid_grid    = total_density > 1e-6;
    exc_fraction  = NaN(n_grid, 1);
    exc_fraction(valid_grid) = exc_density(valid_grid) ./ total_density(valid_grid);

    exc_frac_mat = reshape(exc_fraction, size(GX));

    % Summary metrics
    valid_vals = exc_fraction(valid_grid);
    nw.SpatialEIVariability = std(valid_vals, 'omitnan');

    % Moran's I on grid values (inverse-distance weights between valid grid points)
    if sum(valid_grid) >= 4
        vxy = grid_pts(valid_grid, :);
        D_grid = squareform(pdist(vxy));
        D_grid(1:size(D_grid,1)+1:end) = Inf;
        W_grid = 1 ./ D_grid;
        W_grid(isinf(W_grid)) = 0;
        nw.SpatialEIMoransI = moransI(valid_vals, W_grid);
    end

    network_table = struct2table(nw);

    viz_data.grid_x       = GX;
    viz_data.grid_y       = GY;
    viz_data.exc_fraction = exc_frac_mat;
    viz_data.unit_xy      = xy;
    viz_data.unit_labels  = labels;
end


function I = moransI(x, W)
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
    I = (n / S0) * (z' * (Wv * z)) / (z' * z);
end

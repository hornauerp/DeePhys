function [L_max, cluster_scale, L_curve, r_vals, ci_envelope] = computeRipleysK(xy, study_area, params)
% COMPUTERIPLEYK  Ripley's K/L point pattern analysis.
%
%   [L_max, cluster_scale, L_curve, r_vals, ci_envelope] = computeRipleysK(xy, study_area, params)
%
% Computes the L-function L(r) = sqrt(K(r)/pi) - r, where K(r) counts the
% expected number of additional points within distance r of a random point.
% L(r) > 0 indicates clustering; L(r) < 0 indicates dispersion.
%
% INPUTS:
%   xy          - (N x 2) point coordinates
%   study_area  - scalar area of study region (convex hull area, um^2)
%   params      - struct with fields:
%                   .RipleyMaxRadius   (default 500 um)
%                   .RipleyNBins       (default 50)
%                   .RipleyNSurrogates (default 100) — set 0 to skip CI
%
% OUTPUTS:
%   L_max         - max deviation of L(r) from 0 across all r (summary scalar)
%   cluster_scale - r at which L_max occurs (um)
%   L_curve       - (1 x NRipleyNBins) L(r) values
%   r_vals        - (1 x NRipleyNBins) distance values evaluated
%   ci_envelope   - (2 x NRipleyNBins) [lower; upper] CSR confidence interval

    arguments
        xy          (:,2) double
        study_area  (1,1) double
        params      struct = struct()
    end

    p = mergeDefaults(params, struct( ...
        'RipleyMaxRadius',   500, ...
        'RipleyNBins',       50, ...
        'RipleyNSurrogates', 100));

    n = size(xy, 1);
    r_vals = linspace(0, p.RipleyMaxRadius, p.RipleyNBins + 1);
    r_vals = r_vals(2:end);  % exclude r=0

    L_curve = computeL(xy, r_vals, study_area, n);

    [L_max, idx] = max(abs(L_curve));
    cluster_scale = r_vals(idx);

    ci_envelope = NaN(2, numel(r_vals));
    if p.RipleyNSurrogates > 0 && study_area > 0 && ~isnan(study_area)
        surrogate_L = NaN(p.RipleyNSurrogates, numel(r_vals));
        % Bounding box of convex hull for uniform sampling. Note: using the
        % bounding box while study_area is the convex hull area can inflate
        % L(r) for non-rectangular layouts (surrogates span a larger region).
        x_range = [min(xy(:,1)), max(xy(:,1))];
        y_range = [min(xy(:,2)), max(xy(:,2))];
        for s = 1:p.RipleyNSurrogates
            rand_xy = [rand(n,1) * diff(x_range) + x_range(1), ...
                       rand(n,1) * diff(y_range) + y_range(1)];
            surrogate_L(s,:) = computeL(rand_xy, r_vals, study_area, n);
        end
        ci_envelope(1,:) = prctile(surrogate_L, 2.5, 1);
        ci_envelope(2,:) = prctile(surrogate_L, 97.5, 1);
    end
end


function L = computeL(xy, r_vals, study_area, n)
    D = squareform(pdist(xy));
    D(1:n+1:end) = Inf;  % diagonal to Inf to exclude self
    L = zeros(1, numel(r_vals));
    for ri = 1:numel(r_vals)
        K_r = (study_area / n^2) * sum(D(:) <= r_vals(ri));
        L(ri) = sqrt(K_r / pi) - r_vals(ri);
    end
end


function p = mergeDefaults(p, defaults)
    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        if ~isfield(p, fields{i})
            p.(fields{i}) = defaults.(fields{i});
        end
    end
end

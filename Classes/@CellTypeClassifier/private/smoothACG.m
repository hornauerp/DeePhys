function acg = smoothACG(acg, sigma_bins)
% SMOOTHACG  Gaussian smoothing of ACGs along the lag axis (columns = units).
%
% The centre (zero-lag) bin is excluded from the smoothing input and kept as is: it holds
% no distinct spike pairs and would otherwise bleed into the refractory bins.
    if isempty(acg) || sigma_bins <= 0
        return
    end
    nb = size(acg, 1); c = (nb + 1) / 2;
    h = ceil(4 * sigma_bins);
    k = exp(-0.5 * ((-h:h)' / sigma_bins).^2); k = k / sum(k);
    centre = acg(c, :);
    x = acg; x(c, :) = 0;
    w = ones(nb, 1); w(c) = 0;                  % renormalise where the centre was removed
    s  = conv2(x, k, 'same');
    sw = conv2(w, k, 'same');
    acg = s ./ max(sw, eps);
    acg(c, :) = centre;
end

function response = bootstrapResponse(culture, pre_cutout, post_cutout, bin_size, n_iter, alpha)
% BOOTSTRAPRESPONSE  Bootstrap test for significant firing rate changes between two time windows.
%
% Wraps bootstrapFiringRateResponse using this culture's binned spike matrices.
% Pre and post cutouts must span the same duration.
%
% INPUTS:
%   culture     - Culture object
%   pre_cutout  - [start, end] in seconds — baseline window
%   post_cutout - [start, end] in seconds — post-treatment window
%   bin_size    - Bin width in seconds (default 20)
%   n_iter      - Number of bootstrap iterations (default 1000)
%   alpha       - Two-tailed significance level (default 1e-10)
%
% OUTPUTS:
%   response - struct with fields .increase, .decrease, .unchanged (unit indices)

arguments
    culture     Culture
    pre_cutout  (1,2) double
    post_cutout (1,2) double
    bin_size    (1,1) double = 20
    n_iter      (1,1) double = 1000
    alpha       (1,1) double = 1e-10
end

assert(abs(diff(pre_cutout) - diff(post_cutout)) < 1e-9, ...
    'Pre and post cutouts must span equal durations');

pre_mat  = culture.getBinnedSpikeMat(bin_size, pre_cutout);
post_mat = culture.getBinnedSpikeMat(bin_size, post_cutout);
response = bootstrapFiringRateResponse(pre_mat, post_mat, n_iter, alpha);
end

function [sorted_acgs, sorted_idx] = sortACGsByPeak(acgs, type)
% SORTACGSBYPEAK  Sort an ACG matrix by the peak location of each row.
%
% Sorts autocorrelogram rows by the dominant peak in the first half of each ACG,
% enabling visual grouping by ACG shape (e.g. narrow = fast-spiking interneurons).
%
% INPUTS:
%   acgs - (N_units x N_bins) normalised ACG matrix
%   type - "max"  — sort by bin of maximum amplitude (default)
%           "mean" — sort by mean-weighted peak position
%
% OUTPUTS:
%   sorted_acgs - (N_units x N_bins) reordered ACG matrix
%   sorted_idx  - permutation index used for sorting (for reordering other arrays)

arguments
    acgs double
    type string = "max"
end

half = floor(size(acgs, 2) / 2);

if type == "max"
    [~, sort_key] = max(acgs(:, 1:half), [], 2);
elseif type == "mean"
    sort_key = mean(acgs(:, 1:half), 2);
    sort_key = sort_key ./ max(acgs(:, 1:half), [], 2);
else
    error('Unknown sort type "%s". Use "max" or "mean".', type);
end

[~, sorted_idx] = sort(sort_key, 'descend');
sorted_acgs = acgs(sorted_idx, :);
end

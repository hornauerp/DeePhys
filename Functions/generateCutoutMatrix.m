function cutout_mat = generateCutoutMatrix(rec_length, cutout_idx, cutout)
% GENERATECUTOUTMATRIX  Build absolute time windows for selected recordings in a concatenated timeline.
%
% Given a recording duration and a set of recording indices, computes absolute
% [start, end] times assuming recordings are concatenated back-to-back.
% Useful for extracting specific recordings (e.g. baseline or post-treatment
% time windows) from a concatenated spike matrix.
%
% INPUTS:
%   rec_length  - Duration of each recording in seconds (assumed equal across recordings)
%   cutout_idx  - (1 x N) indices of the recordings to include (1-based)
%   cutout      - (optional) [start, end] offset within each recording in seconds;
%                 default [0, rec_length]
%
% OUTPUTS:
%   cutout_mat  - (N x 2) matrix of absolute [start, end] times in the concatenated timeline

arguments
    rec_length (1,1) double
    cutout_idx (1,:) double
    cutout     (1,2) double = [0, rec_length]
end

cutout_mat = nan(length(cutout_idx), 2);
for i = 1:length(cutout_idx)
    t_start          = cutout(1) + rec_length * (cutout_idx(i) - 1);
    t_end            = cutout(2) + rec_length * (cutout_idx(i) - 1);
    cutout_mat(i,:)  = [t_start, t_end];
end
end

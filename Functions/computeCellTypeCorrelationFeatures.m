function [network_table, unit_table] = computeCellTypeCorrelationFeatures( ...
    sttc_wu, cell_type_labels)
% COMPUTECELLTYPECORRELATIONFEATURES  E/I correlation structure from STTC.
%
%   [network_table, unit_table] = computeCellTypeCorrelationFeatures(sttc_wu, cell_type_labels)
%
% Uses the weighted undirected STTC matrix to compute mean pairwise synchrony
% within and across cell types, and a synchrony index capturing whether
% within-type synchrony exceeds cross-type synchrony.
%
% INPUTS:
%   sttc_wu          - (N x N) symmetric weighted STTC matrix
%   cell_type_labels - (1 x N) 1=excitatory, 2=inhibitory, NaN=unclassified
%
% OUTPUTS:
%   network_table - (1 x 4) table: MeanSTTC_EE, MeanSTTC_II, MeanSTTC_EI, STTCSynchronyIndex
%   unit_table    - (N x 2) table: MeanSTTC_WithinType, MeanSTTC_CrossType

    n_units = numel(cell_type_labels);
    exc_idx = find(cell_type_labels == 1);
    inh_idx = find(cell_type_labels == 2);

    % Network-level: mean pairwise STTC per sub-population
    MeanSTTC_EE = pairMean(sttc_wu, exc_idx, exc_idx, true);
    MeanSTTC_II = pairMean(sttc_wu, inh_idx, inh_idx, true);
    MeanSTTC_EI = pairMean(sttc_wu, exc_idx, inh_idx, false);

    denom = abs(MeanSTTC_EE) + abs(MeanSTTC_EI) + eps;
    STTCSynchronyIndex = (MeanSTTC_EE - MeanSTTC_EI) / denom;

    network_table = table(MeanSTTC_EE, MeanSTTC_II, MeanSTTC_EI, STTCSynchronyIndex);

    % Per-unit: mean STTC with same-type and other-type units
    within_type = NaN(n_units, 1);
    cross_type  = NaN(n_units, 1);

    for u = 1:n_units
        lbl = cell_type_labels(u);
        if isnan(lbl), continue; end
        if lbl == 1
            same_idx  = exc_idx(exc_idx ~= u);
            other_idx = inh_idx;
        else
            same_idx  = inh_idx(inh_idx ~= u);
            other_idx = exc_idx;
        end
        if ~isempty(same_idx)
            within_type(u) = mean(sttc_wu(u, same_idx), 'omitnan');
        end
        if ~isempty(other_idx)
            cross_type(u) = mean(sttc_wu(u, other_idx), 'omitnan');
        end
    end

    unit_table = table(within_type, cross_type, ...
        'VariableNames', {'MeanSTTC_WithinType', 'MeanSTTC_CrossType'});
end


function m = pairMean(W, row_idx, col_idx, upper_tri_only)
% Mean of entries W(row_idx, col_idx), optionally upper-triangle only (for symmetric pairs).
    if isempty(row_idx) || isempty(col_idx)
        m = NaN;
        return
    end
    sub = W(row_idx, col_idx);
    if upper_tri_only && isequal(row_idx, col_idx)
        mask = triu(true(size(sub)), 1);
        vals = sub(mask);
    else
        vals = sub(:);
    end
    m = mean(vals, 'omitnan');
end

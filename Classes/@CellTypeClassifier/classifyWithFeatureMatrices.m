function [labels, train_reduction, test_reduction] = classifyWithFeatureMatrices(ctc, X_train, y_train, X_test, feature_names)
% CLASSIFYWITHFEATUREMATRICES  Supervised UMAP classification from raw feature matrices.
%
% Low-level entry point for mixed or fully external classification workflows.
% Applies joint z-score normalisation (train + test stacked), drops constant
% features, scales by the training-set column maximum, then runs supervisedUMAP.
%
% Normalisation note:
%   DeePhys' default classifyUnits path uses chip-level (ChipID) z-scoring
%   followed by a global z-score. This method uses a single joint z-score —
%   the only option when training and test sets come from different sources.
%
% INPUTS:
%   ctc           - CellTypeClassifier (must have TrainLabels set)
%   X_train       - (N_train x N_features) raw feature matrix for training
%   y_train       - (1 x N_train) class labels (1=excitatory, 2=inhibitory)
%   X_test        - (N_test  x N_features) raw feature matrix for test set
%   feature_names - (1 x N_features) column names from buildFeatureMatrix
%
% OUTPUTS:
%   labels          - (1 x N_test) predicted labels (1=excitatory, 2=inhibitory)
%   train_reduction - (N_train x NDims) UMAP embedding of training set
%   test_reduction  - (N_test  x NDims) UMAP embedding of test set

arguments
    ctc           CellTypeClassifier
    X_train       (:,:) double
    y_train       (1,:) double
    X_test        (:,:) double
    feature_names (1,:) string
end

N_train = size(X_train, 1);

% Joint z-score normalisation
X_combined        = [X_train; X_test];
[X_norm, ~, ~]    = normalize(X_combined, 1, 'zscore');

% Drop features where training block is NaN (constant or all-zero features)
nan_cols = any(isnan(X_norm(1:N_train, :)), 1);
X_norm(:, nan_cols)   = [];
feature_names(nan_cols) = [];

% Scale by training-set column maximum
scale        = max(abs(X_norm(1:N_train, :)), [], 1);
scale(scale == 0) = 1;
X_norm       = X_norm ./ scale;

X_tr = X_norm(1:N_train, :);
X_te = X_norm(N_train+1:end, :);

[labels, ~, ~, test_reduction, train_reduction] = supervisedUMAP(ctc, ...
    X_tr, y_train, feature_names, X_te);
end

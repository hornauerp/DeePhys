function [labels, train_reduction, test_reduction] = classifyWithFeatureMatrices(ctc, X_train, y_train, X_test, feature_names)
% CLASSIFYWITHFEATUREMATRICES  Supervised UMAP classification from raw feature matrices.
%
% Low-level entry point for mixed or fully external classification workflows.
% Applies the shared normalization pipeline (global z-score fitted on training
% data, bad-column removal, training-max scaling) then runs supervisedUMAP.
%
% Per-group z-score (step 1) is the caller's responsibility when required —
% this function handles only steps 2–4 of the normalization pipeline.
%
% INPUTS:
%   ctc           - CellTypeClassifier (must have Parameters set)
%   X_train       - (N_train x N_features) feature matrix for training
%   y_train       - (1 x N_train) class labels (1=excitatory, 2=inhibitory)
%   X_test        - (N_test  x N_features) feature matrix for test set
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

[X_tr, X_te, feature_names] = CellTypeClassifier.normalizeForClassification( ...
    X_train, X_test, feature_names);

[labels, ~, ~, test_reduction, train_reduction] = supervisedUMAP(ctc, ...
    X_tr, y_train, feature_names, X_te);
end

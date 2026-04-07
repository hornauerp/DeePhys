function [Y_pred, umap_train, umap_test, test_reduction, train_reduction, extras] = supervisedUMAP(ctc, X_train, Y_train, feature_names, X_test)
% SUPERVISEDUMAP  Train a supervised UMAP and project test data for classification.
%
% Saves a UMAP template file to ctc.Parameters.UMAP.TemplateDir so that test
% data can be projected onto the same embedding as training data.
% Requires run_umap (MATLAB UMAP toolbox) to be on the MATLAB path.
%
% INPUTS:
%   ctc           - CellTypeClassifier (provides Parameters.UMAP.*)
%   X_train       - (N_train x F) normalised feature matrix
%   Y_train       - (1 x N_train) class labels (1 = excitatory, 2 = inhibitory)
%   feature_names - (1 x F) feature name strings  [kept for API compatibility]
%   X_test        - (N_test x F) normalised feature matrix
%
% OUTPUTS:
%   Y_pred          - (1 x N_test) predicted class labels
%   umap_train      - trained UMAP model object
%   umap_test       - UMAP model used for test projection
%   test_reduction  - (N_test  x NDims) embedding coordinates for test units
%   train_reduction - (N_train x NDims) embedding coordinates for training units
%   extras          - UMAP_extra_results with QF dissimilarity scores (from test call)

arguments
    ctc CellTypeClassifier
    X_train (:,:) double
    Y_train (1,:) double
    feature_names (1,:) string %#ok<INUSA>
    X_test  (:,:) double
end

p = ctc.Parameters.UMAP;
rng(ctc.Parameters.RNGSeed, 'twister');
tdir = p.TemplateDir;
uid           = char(java.util.UUID.randomUUID);
template_file = fullfile(tdir, ['ctc_supervised_umap_template', uid,'.mat']);
clean_up      = onCleanup(@() deleteFile(template_file)); %#ok<NASGU>

% Training call — (label column appended as last column)
train_args = { ...
    'label_column',       'end', ...
    'n_components',       p.SupervisedNDims, ...
    'n_neighbors',        p.SupervisedNNeighbors, ...
    'min_dist',           p.MinDist, ...
    'spread',             p.Spread, ...
    'sgd_tasks',          20, ...
    'method',            'java', ...
    'verbose',            'none', ...
    'save_template_file', template_file};
if isfield(p, 'TargetWeight')
    train_args = [train_args, {'target_weight', p.TargetWeight}];
end
[train_reduction, umap_train, ~, ~] = run_umap([X_train, Y_train'], train_args{:});

% Test call — match_supervisors=3 for label propagation, QF dissimilarity enabled
% Append a dummy label column: the template was trained with label_column='end'
% so it expects F+1 columns. The dummy labels are ignored by match_supervisors.
X_test_with_label = [X_test, zeros(size(X_test, 1), 1)];
test_args = { ...
    'label_column',      'end', ...
    'method',            'java', ...
    'sgd_tasks',         20, ...
    'verbose',           'none', ...
    'cluster_detail',    'adaptive', ...
    'match_supervisors', 3, ...
    'qf_dissimilarity',  true, ...
    'match_scenarios',   4, ...
    'template_file',     template_file};
if isfield(ctc.Parameters, 'Classification') && ...
        isfield(ctc.Parameters.Classification, 'UseJoinedTransform') && ...
        ctc.Parameters.Classification.UseJoinedTransform
    test_args = [test_args, {'joined_transform', true}];
end
[test_reduction, umap_test, ~, extras] = run_umap(X_test_with_label, test_args{:});

Y_pred = extras.supervisorMatchedLabels;
end

function deleteFile(f)
if isfile(f); delete(f); end
end

function [Y_pred, umap_train, umap_test, test_reduction] = supervisedUMAP(ctc, X_train, Y_train, feature_names, X_test)
% SUPERVISEDUMAP  Train a supervised UMAP and project test data for classification.
%
% Saves UMAP template files to ctc.Parameters.UMAP.TemplateDir so that test
% data can be projected onto the same embedding as training data.
% Requires run_umap (MATLAB UMAP toolbox) to be on the MATLAB path.
%
% INPUTS:
%   ctc           - CellTypeClassifier (provides Parameters.UMAP.*)
%   X_train       - (N_train x F) normalised feature matrix
%   Y_train       - (1 x N_train) class labels (1 = excitatory, 2 = inhibitory)
%   feature_names - (1 x F) feature name strings
%   X_test        - (N_test x F) normalised feature matrix
%
% OUTPUTS:
%   Y_pred         - (N_test x 1) predicted class labels
%   umap_train     - trained UMAP model object
%   umap_test      - UMAP model used for test projection
%   test_reduction - (N_test x NDims) embedding coordinates for test units

arguments
    ctc CellTypeClassifier
    X_train (:,:) double
    Y_train (1,:) double
    feature_names (1,:) string
    X_test  (:,:) double
end

p = ctc.Parameters.UMAP;
tdir = p.TemplateDir;
template_file = fullfile(tdir, 'ctc_supervised_umap_template.mat');
train_csv     = fullfile(tdir, 'ctc_umap_train.csv');
test_csv      = fullfile(tdir, 'ctc_umap_test.csv');

% Write training data to CSV (run_umap reads CSV for supervised mode)
train_tbl = array2table([X_train, Y_train'], ...
    'VariableNames', [feature_names, "NeuronType"]);
writetable(train_tbl, train_csv);

[~, umap_train, ~, ~] = run_umap(train_csv, ...
    'label_column',       'end', ...
    'n_components',       p.NDims, ...
    'n_neighbors',        p.NNeighbors, ...
    'min_dist',           p.MinDist, ...
    'spread',             p.Spread, ...
    'sgd_tasks',          20, ...
    'verbose',            'none', ...
    'save_template_file', template_file);

% Write test data to CSV and project using saved template
test_tbl = array2table(X_test, 'VariableNames', feature_names);
writetable(test_tbl, test_csv);

[test_reduction, umap_test, ~, extras] = run_umap(test_csv, ...
    'sgd_tasks',        20, ...
    'verbose',          'none', ...
    'match_supervisors', 3, ...
    'template_file',    template_file);

Y_pred = extras.supervisorMatchedLabels;
end

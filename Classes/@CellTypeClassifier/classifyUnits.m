function classifyUnits(ctc)
% CLASSIFYUNITS  Classify all units using supervised UMAP projection.
%
% Trains a supervised UMAP on the labelled training set, then projects all
% remaining (test) units into the same embedding space to predict cell type.
%
% Requires ctc.TrainLabels to be set (run generateTrainLabels first).
% Sets: ctc.UnitLabels  (1 = excitatory, 2 = inhibitory, NaN = unclassified)

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.TrainLabels), ...
    'Run generateTrainLabels() before classifyUnits()');

p_umap = ctc.Parameters.UMAP;
rg     = ctc.RecordingGroup;
labels = ctc.TrainLabels;

[input_table, ~] = aggregateCultureFeatureTables(rg, "Unit", ...
    p_umap.GroupingVar, p_umap.GroupingValues, 0, [], p_umap.UnitFeatures, [], false);

[X_train, X_test, feature_names] = prepareInputMatrix(rg, input_table, ctc.UnitList, ...
    p_umap.NormalizationVar, labels.umap_train_idx, labels.umap_test_idx);

[Y_pred, ~, ~, test_reduction, train_reduction] = supervisedUMAP(ctc, X_train, labels.sorted_y_train, feature_names, X_test);

ctc.Reduction.Train = train_reduction;
ctc.Reduction.Test  = test_reduction;

full_labels = nan(1, length(ctc.UnitList));
full_labels(labels.sorted_train_ids) = labels.sorted_y_train;
full_labels(labels.umap_test_idx)    = Y_pred;
ctc.UnitLabels = full_labels;

n_exc = sum(full_labels == 1, 'omitnan');
n_inh = sum(full_labels == 2, 'omitnan');
fprintf('Classified %i excitatory, %i inhibitory units (%i unclassified)\n', ...
    n_exc, n_inh, sum(isnan(full_labels)));
end

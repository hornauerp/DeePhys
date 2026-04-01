%% Tutorial 3 — Phenotype Analysis
%
% Covers classification, dimensionality reduction, and regression at unit,
% recording, and culture levels via the Experiment orchestration layer and
% the standalone Classifier / DimReducer / Regressor classes.
%
% Prerequisites:
%   - Tutorial 1 completed (FeatureStore and RecordingProcessors saved)
%   - DeePhys on the MATLAB path

%% 1  Build Experiment from saved processors

proc_paths = {'/path/to/proc1.mat', '/path/to/proc2.mat'};
procs = RecordingProcessor.loadMany(proc_paths);

exp = Experiment.fromProcessors(procs);

% Or from a saved FeatureStore directly (skip the processors):
%   exp = Experiment.fromLegacyGroup(rg);   % from old RecordingGroup

%% 2  Unit-level classification
%
% Classifies each unit using cross-validation grouped by recording,
% so units from the same recording are never split across train/test.

opts = struct();
opts.Algorithm     = 'rf';      % 'rf' (random forest) or 'svm'
opts.KFold         = 5;
opts.FeatureGroups = 'all';     % or ["ActivityFeatures","WaveformFeatures"]
opts.CVLevel       = 'recording';  % group CV at recording level

result_unit = exp.classify('Unit', 'Mutation', opts);

% Inspect result
disp(result_unit);
fprintf('Unit-level accuracy: %.2f ± %.2f\n', ...
    mean(result_unit.Accuracy), std(result_unit.Accuracy));

%% 3  Unit-level classification with parent ACGs
%
% Parent ACGs (from concatenated recording) improve classification by
% providing a fuller spike history. Pass ParentFeatures to use them.

opts_parent = struct();
opts_parent.Algorithm     = 'rf';
opts_parent.FeatureGroups = 'all';
opts_parent.ParentFeatures = 'ACG';   % prefer Parent_ACG* from FeatureStore

result_parent = exp.classify('Unit', 'Mutation', opts_parent);

%% 4  Unit-level dimensionality reduction

opts_umap = struct();
opts_umap.Method       = 'UMAP';
opts_umap.FeatureGroups = 'all';

exp.reduce('Unit', opts_umap);

reduction = exp.Results.DimReduction.Unit.UMAP;
disp(reduction);

% Scatter plot coloured by Mutation
embedding = reduction.Embedding;
labels    = string(exp.FeatureStore.UnitTable.Mutation);
if ~isempty(embedding) && size(embedding, 2) >= 2
    figure;
    gscatter(embedding(:,1), embedding(:,2), labels);
    xlabel('UMAP 1'); ylabel('UMAP 2');
    title('Unit UMAP — coloured by Mutation');
end

%% 5  Recording-level classification

opts_rec = struct();
opts_rec.Algorithm = 'rf';
opts_rec.KFold     = 5;

result_rec = exp.classify('Recording', 'Mutation', opts_rec);
fprintf('Recording-level accuracy: %.2f ± %.2f\n', ...
    mean(result_rec.Accuracy), std(result_rec.Accuracy));

%% 6  Recording-level dimensionality reduction

exp.reduce('Recording', struct('Method', 'PCA'));
pca_result = exp.Results.DimReduction.Recording.PCA;

%% 7  Culture-level classification
%
% Aggregates recordings into culture-level feature vectors (one row per culture).
% Each culture is represented by features at the specified DIV time points,
% concatenated into a wide vector.

opts_cult = struct();
opts_cult.IdentityKeys   = ["ChipID", "PlatingDate"];
opts_cult.GroupingVar    = 'DIV';
opts_cult.GroupingValues = [7, 14, 21, 28];
opts_cult.Algorithm      = 'rf';
opts_cult.KFold          = 5;

result_cult = exp.classify('Culture', 'Mutation', opts_cult);
fprintf('Culture-level accuracy: %.2f ± %.2f\n', ...
    mean(result_cult.Accuracy), std(result_cult.Accuracy));

%% 8  Culture-level regression
%
% Regress a numeric metadata variable (e.g. drug concentration) from
% culture-level features.

opts_reg = struct();
opts_reg.IdentityKeys   = ["ChipID", "PlatingDate"];
opts_reg.GroupingVar    = 'DIV';
opts_reg.GroupingValues = [7, 14, 21, 28];
opts_reg.KFold          = 5;

result_conc = exp.regress('Culture', 'Concentration', opts_reg);
fprintf('Concentration regression R2: %.2f\n', mean(result_conc.R2));

%% 9  Access stored results

% All classification results
disp(fieldnames(exp.Results.Classification));
disp(fieldnames(exp.Results.DimReduction));

% Classification result fields
r = exp.Results.Classification.Mutation;
disp(r);

%% 10  Direct API — bypass Experiment

% Extract feature matrix manually
[X, unit_ids] = exp.FeatureStore.unitMatrix('all');
Y = exp.FeatureStore.UnitTable.Mutation;

% Hierarchy-aware CV: units from the same recording stay together
rec_ids  = exp.FeatureStore.UnitTable.RecordingID;
cv_folds = GroupedCV.byGroups(rec_ids, 5, Y);

% Standalone classification
opts_direct = struct('Algorithm', 'rf', 'KFold', 5);
opts_direct.CVGroups = rec_ids;
result_direct = Classifier.classify(X, Y, opts_direct);

%% 11  Direct API — dimensionality reduction

[X_unit, ~] = exp.FeatureStore.unitMatrix('all');
opts_dr = struct('Method', 'UMAP', 'NDims', 2);
dim_result = DimReducer.reduce(X_unit, opts_dr);
disp(dim_result);

%% 12  Within-group normalization
%
% Normalize features per recording (or per ChipID) before ML to remove
% batch effects. Specify the metadata field to group by.

opts_norm = struct();
opts_norm.NormalizationVar = 'ChipID';
opts_norm.Algorithm        = 'rf';

result_norm = exp.classify('Unit', 'Mutation', opts_norm);

%% 13  Feature group selection for classification

opts_fg = struct();
opts_fg.FeatureGroups = ["ActivityFeatures", "WaveformFeatures"];
opts_fg.Algorithm     = 'rf';

result_fg = exp.classify('Unit', 'Mutation', opts_fg);
fprintf('Activity+Waveform accuracy: %.2f\n', mean(result_fg.Accuracy));

%% 14  Confusion matrix visualization

if ~isempty(result_unit.ConfusionMatrix)
    figure;
    confusionchart(result_unit.ConfusionMatrix, result_unit.ClassNames);
    title('Unit classification — confusion matrix');
end

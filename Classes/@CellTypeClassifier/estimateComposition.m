function T = estimateComposition(ctc, opts)
% ESTIMATECOMPOSITION  Inhibitory fraction per culture and identity per unit from calibrated
% probabilities (quantification route; complements the graph labels of classifyUnits).
%
%   T = ctc.estimateComposition()
%   T = ctc.estimateComposition(GroupBy="ChipID", NComponents=30, NBoot=200)
%
% Uses the same harmonised, normalised features and the same training seeds as classifyUnits
% (run generateTrainLabels first). A ridge logistic regression on the leading principal
% components gives each unit a probability of being inhibitory (seeds: out-of-fold). Each
% culture's inhibitory fraction is then re-estimated by EM under prior shift and the unit
% probabilities are adjusted to it (estimateClassComposition / adjustClassPrior), so counting
% labels no longer pulls every culture towards the 50 % of the balanced seeds.
%
% OPTIONS:
%   GroupBy      - UnitTable field(s) defining a culture (default Parameters.CultureKeys)
%   NComponents  - principal components of the feature matrix (default 30)
%   NFolds, NBoot, Lambda - see estimateClassComposition
%   Training     - "auto" (default): every labelled unit when both classes are labelled
%                  (GroundTruthMethod "metadata", e.g. pure cultures), else the graph seeds;
%                  "seeds" | "labelled" to force one. The graph seeds keep only the
%                  most typical subset of the labelled units.
%
% OUTPUT / SETS:
%   T                         - table per culture: Group, NUnits, Fraction, CILow, CIHigh,
%                               Counted (share p > 0.5), MeanProbability
%   ctc.Composition.Table     - T
%   ctc.Composition.Probability (1 x N) adjusted probability of being inhibitory per unit row
%   ctc.Composition.RawProbability (1 x N) probability at the seed prior
%   ctc.Composition.Labels    (1 x N) 1 = excitatory, 2 = inhibitory (Probability > 0.5)
%   ctc.Composition.NSpikes   (1 x N) spikes behind each unit's features
%   ctc.Composition.Model     PCA (Coeff, Mu) and logistic model, to score new feature rows
%                             (spikeCountReliability)
% Keep only confident units for downstream analyses with ctc.confidentLabels().
arguments
    ctc CellTypeClassifier
    opts.GroupBy string = string.empty
    opts.NComponents (1,1) double = 30
    opts.NFolds (1,1) double = 5
    opts.NBoot  (1,1) double = 200
    opts.Lambda double = []
    opts.Training (1,1) string {mustBeMember(opts.Training, ["auto", "seeds", "labelled"])} = "auto"
end
assert(~isempty(ctc.TrainLabels) && ~isempty(ctc.NormalizationParams), ...
    'CellTypeClassifier:noTrainLabels', 'Run generateTrainLabels() before estimateComposition().');

% -- features exactly as in classifyUnits ---------------------------------------------
ctc.buildNormalizedFeatures();
nf = ctc.NormalizedFeatures;
np = ctc.NormalizationParams;
X = normalize(nf.X_pergroup, 'center', np.mu_global, 'scale', np.sigma_global);
X(:, np.nan_cols) = [];
X = X ./ np.scale;
if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X = X(:, np.feature_selection_mask);
end
X(~isfinite(X)) = 0;
n_pc = min([opts.NComponents, size(X, 2), size(X, 1) - 1]);
[coeff, Z, ~, ~, ~, mu_pca] = pca(X, 'NumComponents', n_pc);

% -- cultures of the unique units ------------------------------------------------------
ut = ctc.FeatureStore.UnitTable;
if isempty(opts.GroupBy)
    cult = FeatureStore.getCultureIDsForUnits(ut.RecordingID, ctc.FeatureStore.MetadataTable, ...
        ctc.Parameters.CultureKeys);
else
    cult = join(string(ut{:, cellstr(opts.GroupBy)}), "|", 2);
end
cult = string(cult(:));
groups = cult(nf.unique_to_fs);

% -- training units: the graph seeds, or every labelled unit when both classes carry explicit
% labels (pure-culture metadata). Pure cultures label all their units; the enrichment seeds
% keep only the most typical subset of them, which can bias the boundary.
labels = ctc.TrainLabels;
train_ids = labels.sorted_train_ids(:); y_train = labels.sorted_y_train(:);
has_ce = ~isempty(ctc.CounterexampleUnitIdx) && any(ctc.CounterexampleUnitIdx);
use_labelled = opts.Training == "labelled" || (opts.Training == "auto" && has_ce);
if use_labelled
    assert(has_ce, 'CellTypeClassifier:noExplicitLabels', ...
        'Training "labelled" needs explicit labels for both classes (GroundTruthMethod "metadata").');
    resp = ctc.ResponsiveUnitIdx(:); ce = ctc.CounterexampleUnitIdx(:);
    resp_u = resp(nf.unique_to_fs); ce_u = ce(nf.unique_to_fs);
    r_lab = ctc.Parameters.TrainLabels.ResponsiveClassLabel;
    train_ids = find(resp_u | ce_u);
    y_train = r_lab * resp_u(train_ids) + (3 - r_lab) * (ce_u(train_ids) & ~resp_u(train_ids));
end
[T, post, p, mdl] = estimateClassComposition(Z, train_ids, y_train, groups, ...
    NFolds=opts.NFolds, NBoot=opts.NBoot, Lambda=opts.Lambda, Seed=ctc.Parameters.RNGSeed);

a2u = nf.all_to_unique;
training = "seeds"; if use_labelled, training = "labelled"; end
% spikes behind each unit's features (its representative's spike train)
n_spk = arrayfun(@(u) numel(u.SpikeTimes), nf.unique_ud);
ctc.Composition = struct('Table', T, ...
    'Probability', post(a2u)', 'RawProbability', p(a2u)', ...
    'Labels', 1 + double(post(a2u)' > 0.5), ...
    'NSpikes', n_spk(a2u), ...
    'TrainPrior', mean(y_train == 2), 'NComponents', n_pc, ...
    'Training', training, 'NTrain', numel(train_ids), ...
    'Model', struct('Coeff', coeff, 'Mu', mu_pca, 'Classifier', mdl));
fprintf(['estimateComposition: %d cultures, inhibitory fraction median %.1f%% [%.1f %.1f] ' ...
    '(counting p > 0.5 would give %.1f%%)\n'], height(T), 100 * median(T.Fraction), ...
    100 * prctile(T.Fraction, 25), 100 * prctile(T.Fraction, 75), 100 * median(T.Counted));
end

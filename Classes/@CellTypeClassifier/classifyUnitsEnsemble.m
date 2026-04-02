function classifyUnitsEnsemble(ctc)
% CLASSIFYUNITSENSEMBLE  Classify units using majority vote across multiple RNG seeds.
%
% Runs the full generateTrainLabels + classifyUnits pipeline N times with
% different RNG seeds and assigns majority-vote labels. This directly addresses
% the dominant source of non-reproducibility in the pipeline (sensitivity to a
% single seed choice).
%
% Confidence = fraction of seeds agreeing on each unit's label (0-1).
% Units that fail to reach MinAgreement across seeds are set to NaN.
%
% Parameters are read from ctc.Parameters.Ensemble:
%   Enabled      - must be true to call this method (checked internally)
%   Seeds        - (1 x N) RNG seeds, one per run (default [42,1042,2042,3042,4042])
%   MinAgreement - minimum vote fraction for label assignment (default 0.6)
%
% Requires ctc.ResponsiveUnitIdx to be set (run identifyResponsiveUnits first).
%
% Sets:
%   ctc.UnitLabels     - majority-vote labels (1=exc, 2=inh, NaN=below threshold)
%   ctc.UnitConfidence - fraction of seeds agreeing on each unit's label

arguments
    ctc CellTypeClassifier
end

assert(~isempty(ctc.ResponsiveUnitIdx), ...
    'Run identifyResponsiveUnits() before classifyUnitsEnsemble().');

p_ens   = ctc.Parameters.Ensemble;
seeds   = p_ens.Seeds;
N       = numel(seeds);
n_units = numel(ctc.UnitDataArray);

orig_seed = ctc.Parameters.RNGSeed;

label_matrix = nan(N, n_units);

fprintf('Ensemble classification: %d runs (seeds: %s, MinAgreement=%.2f)\n', ...
    N, num2str(seeds), p_ens.MinAgreement);

for r = 1:N
    fprintf('  Run %d/%d (seed=%d)...\n', r, N, seeds(r));
    ctc.Parameters.RNGSeed = seeds(r);

    try
        ctc.generateTrainLabels();
        ctc.classifyUnits();
        label_matrix(r, :) = ctc.UnitLabels;
    catch ME
        warning('CellTypeClassifier:ensembleRunFailed', ...
            'Run %d (seed=%d) failed: %s', r, seeds(r), ME.message);
    end
end

% Restore original seed
ctc.Parameters.RNGSeed = orig_seed;

% -- Majority vote -----------------------------------------------------------
final_labels     = nan(1, n_units);
final_confidence = zeros(1, n_units);

for u = 1:n_units
    votes       = label_matrix(:, u);
    valid_votes = votes(~isnan(votes));
    if isempty(valid_votes)
        continue
    end
    [mode_label, freq] = mode(valid_votes);
    agreement = freq / N;
    if agreement >= p_ens.MinAgreement
        final_labels(u)     = mode_label;
        final_confidence(u) = agreement;
    end
    % Units below MinAgreement remain NaN with confidence = max agreement seen
    if isnan(final_labels(u))
        final_confidence(u) = freq / N;
    end
end

ctc.UnitLabels     = final_labels;
ctc.UnitConfidence = final_confidence;

n_exc = sum(final_labels == 1, 'omitnan');
n_inh = sum(final_labels == 2, 'omitnan');
n_unc = sum(isnan(final_labels));
fprintf('Ensemble result: %d excitatory, %d inhibitory, %d unclassified (agreement < %.2f)\n', ...
    n_exc, n_inh, n_unc, p_ens.MinAgreement);
end

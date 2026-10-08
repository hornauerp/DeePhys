function [labels, info] = confidentLabels(ctc, opts)
% CONFIDENTLABELS  Unit labels restricted to confident predictions (NaN for the others).
%
%   labels = ctc.confidentLabels()                          % P(inhibitory) >= 0.8 or <= 0.2
%   [labels, info] = ctc.confidentLabels(MinProbability=0.9, MinSpikes=300)
%
% For downstream analyses that should rest on confident units only (e.g. EIAnalyzer with
% Labels.Source = "confident"). Two sources of P(inhibitory) per unit:
%   "composition" (default when estimateComposition has run): the calibrated probability
%       given the unit's features and its culture's estimated composition. If calibrated, the
%       expected share of wrong labels among the kept units is mean(1 - max(P, 1 - P)),
%       reported as info.ExpectedError.
%   "graph": UnitLabels / UnitConfidence of classifyUnits (label-propagation class
%       probability). It ranks units but is not calibrated, so no error rate is reported.
% MinSpikes removes units whose features rest on fewer spikes. spikeCountReliability shows how
% stable labels are under spike-train thinning (pass MinSpikes=R.MinSpikes to apply its value;
% no spike filter is applied by default).
%
% OPTIONS:
%   Source         - "auto" (composition if available, else graph) | "composition" | "graph"
%   MinProbability - keep units with P >= MinProbability or P <= 1 - MinProbability (0.8);
%                    0.5 keeps every unit
%   MinSpikes      - minimum spikes behind the features (default 0)
%
% OUTPUT:
%   labels - (1 x N) 1 = excitatory, 2 = inhibitory, NaN = not confident / too few spikes
%   info   - struct: Source, MinProbability, MinSpikes, Kept (1 x N), Probability (1 x N),
%            NSpikes (1 x N), Coverage, NExcitatory, NInhibitory, ExpectedError
arguments
    ctc CellTypeClassifier
    opts.Source (1,1) string {mustBeMember(opts.Source, ["auto", "composition", "graph"])} = "auto"
    opts.MinProbability (1,1) double {mustBeInRange(opts.MinProbability, 0.5, 1)} = 0.8
    opts.MinSpikes double = []
end
has_comp = isfield(ctc.Composition, 'Probability') && ~isempty(ctc.Composition.Probability);
source = opts.Source;
if source == "auto"
    source = "graph"; if has_comp, source = "composition"; end
end
if source == "composition"
    assert(has_comp, 'CellTypeClassifier:noComposition', 'Run estimateComposition() first.');
    p_inh = ctc.Composition.Probability(:)';
    n_spk = ctc.Composition.NSpikes(:)';
else
    assert(~isempty(ctc.UnitLabels), 'CellTypeClassifier:notClassified', 'Run classifyUnits() first.');
    lab = ctc.UnitLabels(:)'; conf = ctc.UnitConfidence(:)';
    p_inh = conf; p_inh(lab == 1) = 1 - conf(lab == 1); p_inh(isnan(lab)) = NaN;
    ctc.buildNormalizedFeatures();
    nf = ctc.NormalizedFeatures;
    n_u = arrayfun(@(u) numel(u.SpikeTimes), nf.unique_ud);
    n_spk = n_u(nf.all_to_unique);
end
min_spk = opts.MinSpikes;
if isempty(min_spk), min_spk = 0; end
[labels, kept, exp_err] = CellTypeClassifier.selectConfident(p_inh, n_spk, opts.MinProbability, min_spk);
if source == "graph", exp_err = NaN; end
info = struct('Source', source, 'MinProbability', opts.MinProbability, 'MinSpikes', min_spk, ...
    'Kept', kept, 'Probability', p_inh, 'NSpikes', n_spk, 'Coverage', mean(kept), ...
    'NExcitatory', sum(labels == 1), 'NInhibitory', sum(labels == 2), 'ExpectedError', exp_err);
fprintf(['confidentLabels (%s, P >= %.2f either way, >= %d spikes): %d of %d units kept (%.0f%%), ' ...
    '%d excitatory / %d inhibitory'], source, opts.MinProbability, min_spk, sum(kept), numel(kept), ...
    100 * mean(kept), info.NExcitatory, info.NInhibitory);
if isfinite(exp_err)
    fprintf('; expected share of wrong labels %.1f%%\n', 100 * exp_err);
else
    fprintf('\n');
end
end

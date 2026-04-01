function identifyResponsiveUnits(ctc, metadata_filter)
% IDENTIFYRESPONSIVEUNITS  Identify units with significant firing rate changes.
%
% Supports two ground truth methods (set via Parameters.Bootstrap.GroundTruthMethod):
%
%   "full_curve" (default) — For each culture, computes mean firing rate per
%       recording (concentration level) and tests for monotonic increase using
%       Spearman rank correlation with a permutation test. Cultures with fewer
%       than MinRecordings recordings fall back to the two-window method.
%       Strength = absolute Spearman rho (0-1).
%
%   "two_window" — Legacy method comparing pre-treatment vs post-treatment
%       firing rates using a bootstrap permutation test.
%       Strength = |empirical difference| / null std (effect size).
%
% Units that significantly increase are flagged as inhibitory candidates
% (ResponsiveUnitIdx = true). For each unit, a continuous strength score is
% stored in ctc.ResponsiveStrength — higher values indicate more credible
% responsive identification.
%
% When MinResponsiveStrength > 0, responsive units below this strength
% threshold are demoted to non-responsive. This filters out units that
% barely pass the significance test but have weak dose-response trends.
%
% INPUTS:
%   ctc             - CellTypeClassifier object
%   metadata_filter - (optional) {field, value} cell pair to restrict which cultures
%                     are tested (e.g. {'AAV', 128}). Default: all cultures.
%
% Sets: ctc.ResponsiveUnitIdx, ctc.ResponsiveStrength, ctc.UnitList

arguments
    ctc CellTypeClassifier
    metadata_filter cell = {}   % e.g. {'AAV', 128}
end

p = ctc.Parameters.Bootstrap;
rg = ctc.RecordingGroup;
method = p.GroundTruthMethod;

responsive_idx = [];
strength_all   = [];
unit_list = Unit.empty;

for c = 1:length(rg.Cultures)
    culture = rg.Cultures(c);

    % Apply optional metadata filter — only count increases from matching cultures
    if ~isempty(metadata_filter) && isfield(culture.Metadata, metadata_filter{1})
        culture_matches = any([culture.Metadata.(metadata_filter{1})] == metadata_filter{2});
    else
        culture_matches = true;
    end

    % Select ground truth method
    if method == "full_curve" && culture.N >= p.MinRecordings
        response = culture.bootstrapDoseResponse(p.NIter, p.Alpha);
    else
        if method == "full_curve" && culture.N < p.MinRecordings
            fprintf('Culture %d has only %d recordings (<%d) — falling back to two_window\n', ...
                c, culture.N, p.MinRecordings);
        end
        response = culture.bootstrapResponse(p.PreCutout, p.PostCutout, p.BinSize, p.NIter, p.Alpha);
    end

    n_prev = length(unit_list);
    if culture_matches
        responsive_idx = [responsive_idx, response.increase + n_prev]; %#ok<AGROW>
    end
    strength_all = [strength_all, response.strength]; %#ok<AGROW>
    unit_list = [unit_list, culture.Units]; %#ok<AGROW>
end

ctc.UnitList = unit_list;
ctc.ResponsiveUnitIdx = false(1, length(unit_list));
ctc.ResponsiveUnitIdx(responsive_idx) = true;
ctc.ResponsiveStrength = strength_all;

% ── Optional strength-based filtering ────────────────────────────────────────
if isfield(p, 'MinResponsiveStrength') && p.MinResponsiveStrength > 0
    weak = ctc.ResponsiveUnitIdx & (strength_all < p.MinResponsiveStrength);
    n_demoted = sum(weak);
    ctc.ResponsiveUnitIdx(weak) = false;
    if n_demoted > 0
        fprintf('Demoted %d weak responders (strength < %.2f)\n', ...
            n_demoted, p.MinResponsiveStrength);
    end
end

fprintf('Inhibitory candidates: %i / %i units (%.1f%%)\n', ...
    sum(ctc.ResponsiveUnitIdx), length(ctc.UnitList), ...
    100 * mean(ctc.ResponsiveUnitIdx));
end

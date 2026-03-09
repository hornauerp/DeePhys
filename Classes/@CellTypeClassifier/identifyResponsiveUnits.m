function identifyResponsiveUnits(ctc, metadata_filter)
% IDENTIFYRESPONSIVEUNITS  Bootstrap test to identify units with significant firing rate increases.
%
% For each culture, compares pre-treatment vs post-treatment firing rates using
% a bootstrap permutation test. Units that significantly increase are flagged as
% inhibitory candidates (ResponsiveUnitIdx = true).
%
% INPUTS:
%   ctc             - CellTypeClassifier object
%   metadata_filter - (optional) {field, value} cell pair to restrict which cultures
%                     are tested (e.g. {'AAV', 128} to only flag responders from
%                     cultures with AAV=128). Default: all cultures.
%
% Sets: ctc.ResponsiveUnitIdx, ctc.UnitList

arguments
    ctc CellTypeClassifier
    metadata_filter cell = {}   % e.g. {'AAV', 128}
end

p = ctc.Parameters.Bootstrap;
rg = ctc.RecordingGroup;

responsive_idx = [];
unit_list = Unit.empty;

for c = 1:length(rg.Cultures)
    culture = rg.Cultures(c);

    % Bootstrap test: bin spike matrix and compare pre vs post windows
    pre_mat  = culture.getBinnedSpikeMat(p.BinSize, p.PreCutout);
    post_mat = culture.getBinnedSpikeMat(p.BinSize, p.PostCutout);
    response = bootstrapFiringRateResponse(pre_mat, post_mat, p.NIter, p.Alpha);

    % Apply optional metadata filter — only count increases from matching cultures
    if ~isempty(metadata_filter) && isfield(culture.Metadata, metadata_filter{1})
        culture_matches = any([culture.Metadata.(metadata_filter{1})] == metadata_filter{2});
    else
        culture_matches = true;
    end

    if culture_matches
        responsive_idx = [responsive_idx, response.increase + length(unit_list)]; %#ok<AGROW>
    end
    unit_list = [unit_list, culture.Units]; %#ok<AGROW>
end

ctc.UnitList = unit_list;
ctc.ResponsiveUnitIdx = false(1, length(unit_list));
ctc.ResponsiveUnitIdx(responsive_idx) = true;
fprintf('Inhibitory candidates: %i / %i units (%.1f%%)\n', ...
    sum(ctc.ResponsiveUnitIdx), length(ctc.UnitList), ...
    100 * mean(ctc.ResponsiveUnitIdx));
end

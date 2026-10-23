classdef FeatureCatalog
% FEATURECATALOG  Central registry of named feature sets.
%
% Defines which features belong to each (group, set) combination.
% Used by FeatureStore to filter columns at retrieval time. All methods
% are Static — no instance needed.
%
% BUILT-IN SETS:
%   "full" — all features in each group (default, backwards-compatible)
%   "core" — curated non-redundant features with clear biological meaning:
%             • Non-redundant: no algebraically determined duplicates
%             • Interpretable: each feature has a clear neuroscientific meaning
%             • Stable: low sensitivity to edge cases (e.g. unstable AUC features removed)
%
% CUSTOM SETS:
%   Users can register their own named sets for the current session:
%     FeatureCatalog.registerSet("pub", ["FiringRate","T2Pdelay","HalfWidth"])
%     FeatureCatalog.removeSet("pub")
%   Custom sets are cleared when the class is cleared (clear FeatureCatalog).
%   Set names must be valid MATLAB identifiers (isvarname).
%   Built-in sets "core" and "full" cannot be overwritten or removed.
%
% USAGE:
%   names = FeatureCatalog.features("ActivityFeatures", "core")
%   names = FeatureCatalog.features("all", "core")   % union across all groups
%   sets  = FeatureCatalog.availableSets()            % → ["core","full",<custom...>]
%   tf    = FeatureCatalog.isKnownFeature("FiringRate")  % → true
%   FeatureCatalog.registerSet("pub", ["FiringRate","T2Pdelay"])
%   FeatureCatalog.removeSet("pub")

    methods (Static)

        function names = features(group, set)
        % FEATURES  Return feature names for a (group, set) combination.
        %
        %   names = FeatureCatalog.features("ActivityFeatures", "core")
        %   names = FeatureCatalog.features("all", "core")
        %   names = FeatureCatalog.features("all", "full")  % → string.empty (no filter)
        %   names = FeatureCatalog.features("all", "pub")   % custom set
        %
        % For "full" set with group="all", returns string.empty as a sentinel
        % meaning "no column filter" — FeatureStore.selectFeatureCols skips
        % filtering when the whitelist is empty.
            arguments
                group (1,1) string
                set   (1,1) string = "full"
            end
            valid = FeatureCatalog.availableSets();
            if ~ismember(set, valid)
                error('FeatureCatalog:unknownSet', ...
                    'Unknown feature set "%s". Valid sets: %s.', set, strjoin(valid, ", "));
            end
            switch set
                case "full"
                    if group == "all"
                        % Sentinel: no filter applied — return empty
                        names = string.empty;
                    else
                        names = FeatureCatalog.fullSet(group);
                    end
                case "core"
                    if group == "all"
                        % Union of core sets across all known groups
                        groups = FeatureCatalog.groups();
                        names = string.empty;
                        for g = groups
                            names = [names, FeatureCatalog.coreSet(g)]; %#ok<AGROW>
                        end
                        names = unique(names, 'stable');
                    else
                        names = FeatureCatalog.coreSet(group);
                    end
                otherwise
                    % Custom set: explicit name list registered via registerSet().
                    % group="all" → return all registered names.
                    % group=<specific> → return intersection with that group's full set.
                    reg       = FeatureCatalog.registry('get');
                    all_names = reg.(char(set));
                    if group == "all"
                        names = all_names;
                    else
                        group_full = FeatureCatalog.fullSet(group);
                        if isempty(group_full)
                            names = string.empty;
                        else
                            names = all_names(ismember(all_names, group_full));
                        end
                    end
            end
        end

        function g = groups()
        % GROUPS  Names of the scalar feature groups (the vector groups ReferenceWaveform /
        %   ACG / FullACG are listed by rawGroups).
            g = ["ActivityFeatures","WaveformFeatures","RegularityFeatures", ...
                 "Catch22","BombcellMetrics","GraphFeatures","BurstParticipation", ...
                 "NetworkBurst","NetworkRegularity","NetworkCatch22","NetworkGraph", ...
                 "CellTypeGraph","CellTypeBalance", ...
                 "CellTypeBurst","CellTypeActivity","CellTypeCorrelation", ...
                 "SpatialFeatures"];
        end

        function g = rawGroups()
        % RAWGROUPS  Vector-valued column groups (Waveform1..N, ACG1..N, FullACG...).
            g = ["ReferenceWaveform", "Waveform", "ACG", "FullACG"];
        end

        function tf = isGroup(name)
        % ISGROUP  True for a known (scalar or vector) feature group name.
            tf = ismember(string(name), [FeatureCatalog.groups(), FeatureCatalog.rawGroups()]);
        end

        function tf = matchColumns(group, cols)
        % MATCHCOLUMNS  Which of the column names cols belong to a feature group: the single
        %   source of group membership for FeatureStore. Exact names, except
        %     - graph groups: stems resolve to their algorithm-suffixed columns
        %       (Density -> Density_CCG, Density_STTC_w; RichClub -> RichClub1_CCG);
        %       CellTypeGraph takes the _EE / _II columns, the other graph groups do not;
        %     - SpatialFeatures: base names and their E/I-stratified variants (_E, _I, _EE, ...);
        %     - vector groups: Waveform1..N, ACG1..N, FullACG*.
        %   An unknown group matches nothing.
            cols = string(cols);
            group = string(group);
            switch group
                case {"ReferenceWaveform", "Waveform"}
                    tf = FeatureCatalog.matchRegexp(cols, '^Waveform\d+$');
                case "ACG"
                    tf = FeatureCatalog.matchRegexp(cols, '^ACG\d+$');
                case "FullACG"
                    tf = startsWith(cols, "FullACG");
                case "Catch22"
                    tf = ismember(cols, FeatureCatalog.fullSet("Catch22"));
                    if ~any(tf), tf = startsWith(cols, "SC_"); end
                case {"GraphFeatures", "NetworkGraph"}
                    tf = FeatureCatalog.matchStems(cols, FeatureCatalog.fullSet(group)) & ...
                         ~endsWith(cols, ["_EE", "_II"]);
                case "CellTypeGraph"
                    roots = [FeatureCatalog.fullSet("NetworkGraph"), FeatureCatalog.fullSet("GraphFeatures")];
                    tf = (FeatureCatalog.matchStems(cols, roots) & endsWith(cols, ["_EE", "_II"])) | ...
                         startsWith(cols, ["CrossTypeDensity", "WithinTypeConnections", "CrossTypeConnections"]);
                case "SpatialFeatures"
                    base = FeatureCatalog.fullSet("SpatialFeatures");
                    variants = base(:) + ["_E", "_I", "_EE", "_II", "_EI"];
                    tf = ismember(cols, base) | ismember(cols, variants(:));
                otherwise
                    if FeatureCatalog.isGroup(group)
                        tf = ismember(cols, FeatureCatalog.fullSet(group));
                    else
                        tf = false(size(cols));
                    end
            end
        end

        function names = columns(group, cols)
        % COLUMNS  The column names among cols that belong to a group (see matchColumns), e.g.
        %   FeatureCatalog.columns("NetworkGraph", X.Properties.VariableNames) for
        %   Classifier.classifyByGroups, which needs exact column names.
            cols = string(cols);
            names = cols(FeatureCatalog.matchColumns(group, cols));
        end

        function tf = inSet(cols, set)
        % INSET  Which columns belong to a named feature set ("core" or a custom set). Names
        %   match exactly; graph stems (Density, Modularity, ...) also match their
        %   algorithm-suffixed columns. "full" keeps every column.
            cols = string(cols);
            if string(set) == "full"
                tf = true(size(cols));
                return
            end
            whitelist = FeatureCatalog.features("all", set);
            graph_stems = [FeatureCatalog.fullSet("GraphFeatures"), FeatureCatalog.fullSet("NetworkGraph"), ...
                           FeatureCatalog.fullSet("CellTypeGraph")];
            tf = ismember(cols, whitelist) | ...
                 FeatureCatalog.matchStems(cols, whitelist(ismember(whitelist, graph_stems)));
        end

        function sets = availableSets()
        % AVAILABLESETS  Return the names of all available feature sets.
        %   Built-in sets appear first; custom sets follow in registration order.
            reg    = FeatureCatalog.registry('get');
            custom = string(fieldnames(reg))';   % (1,:) string, registration order
            sets   = ["core", "full", custom];
        end

        function tf = isKnownFeature(name)
        % ISKNOWNFEATURE  True if name appears in any group's full set.
        %   Useful for distinguishing scalar features from raw vector columns
        %   (Waveform1..N, ACG1..N) which are not in any named group.
            arguments
                name (1,1) string
            end
            tf = false;
            for g = FeatureCatalog.groups()
                if ismember(name, FeatureCatalog.fullSet(g))
                    tf = true;
                    return
                end
            end
        end

        function registerSet(name, feature_names)
        % REGISTERSET  Register a named custom feature set for this session.
        %
        %   FeatureCatalog.registerSet("pub", ["FiringRate","T2Pdelay","HalfWidth"])
        %
        % Rules:
        %   • name must be a valid MATLAB identifier (isvarname).
        %   • Built-in sets "core" and "full" cannot be overwritten.
        %   • Re-registering an existing custom set silently replaces it.
        %   • feature_names may be string.empty (produces a set that returns no columns).
        %   • Duplicate names within feature_names are collapsed (unique, stable order).
        %
        % Custom sets are session-scoped: cleared by "clear FeatureCatalog".
            arguments
                name          (1,1) string
                feature_names (1,:) string = string.empty
            end
            if ~isvarname(name)
                error('FeatureCatalog:invalidName', ...
                    '"%s" is not a valid MATLAB identifier. Set names must satisfy isvarname().', name);
            end
            if ismember(name, ["core","full"])
                error('FeatureCatalog:builtinImmutable', ...
                    'Built-in set "%s" cannot be overwritten.', name);
            end
            FeatureCatalog.registry('set', char(name), unique(feature_names, 'stable'));
        end

        function removeSet(name)
        % REMOVESET  Unregister a previously registered custom feature set.
        %
        %   FeatureCatalog.removeSet("pub")
        %
        % Errors if name is a built-in set or is not currently registered.
            arguments
                name (1,1) string
            end
            if ismember(name, ["core","full"])
                error('FeatureCatalog:builtinImmutable', ...
                    'Built-in set "%s" cannot be removed.', name);
            end
            reg = FeatureCatalog.registry('get');
            if ~isfield(reg, char(name))
                error('FeatureCatalog:unknownSet', ...
                    'No custom set named "%s" is registered. Call availableSets() to see what exists.', name);
            end
            FeatureCatalog.registry('remove', char(name));
        end

    end

    % =========================================================================
    % Private set definitions and session registry
    % =========================================================================
    methods (Static, Access = private)

        function tf = matchStems(cols, stems)
        % Columns named stem, stem<digits> or stem[<digits>]_<suffix>; stems ending in "_"
        % (Strength_, MeanStrength_) match any continuation.
            tf = false(size(cols));
            for s = string(stems(:))'
                if endsWith(s, "_")
                    tf = tf | startsWith(cols, s);
                else
                    tf = tf | FeatureCatalog.matchRegexp(cols, ['^' regexptranslate('escape', char(s)) '\d*(_.*)?$']);
                end
            end
        end

        function tf = matchRegexp(cols, pattern)
            tf = reshape(~cellfun(@isempty, regexp(cellstr(cols), pattern, 'once')), size(cols));
        end

        function reg = registry(action, name, value)
        % REGISTRY  Read/write the session-scoped custom set store.
        %
        %   reg = FeatureCatalog.registry('get')
        %   FeatureCatalog.registry('set',    name, value)
        %   FeatureCatalog.registry('remove', name)
        %   FeatureCatalog.registry('reset')          % clears all (for testing)
        %
        % Backed by a persistent struct; lives until "clear FeatureCatalog".
        % Field names = set names; values = string (1,:) arrays of feature names.
            persistent store
            if isempty(store)
                store = struct();
            end
            switch action
                case 'get'
                    reg = store;
                case 'set'
                    store.(name) = value;
                    reg = store;
                case 'remove'
                    store = rmfield(store, name);
                    reg = store;
                case 'reset'
                    store = struct();
                    reg = store;
                otherwise
                    error('FeatureCatalog:registry', 'Unknown action "%s".', action);
            end
        end

        function names = coreSet(group)
        % CORESET  Non-redundant, interpretable feature names per group.
        %
        % Curation rationale:
        %   ActivityFeatures: MeanISI = 1/FR (redundant), VarianceISI = f(CV,MeanISI)
        %     (redundant), CV is dominated by burst structure and strictly less
        %     informative than CV2/LvR, LV has lower refractory-period correction
        %     than LvR, PartialAutocorrelation at lag-1 captures little beyond CV2.
        %   WaveformFeatures: T2Pdelay is the primary E/I discriminator; HalfWidth
        %     and Asymmetry add complementary shape information. AUC features depend
        %     on unstable zero-crossing detection. Rise/Decay use slewrate() on AP
        %     waveform segments (semantic mismatch). T2Pratio correlates with T2Pdelay.
        %   RegularityFeatures: Magnitude is correlated with firing rate (normalization
        %     artefact); Frequency and Fit capture the independent variance.
        %   Catch22: generic time-series features; overlaps strongly with activity
        %     features; excluded from core to keep it interpretable.
        %   BombcellMetrics: fractionRPVs, waveformDuration, presenceRatio are the
        %     three most diagnostically independent QC metrics.
        %   GraphFeatures: excluded from core — meaningful only for dense networks;
        %     NaN-dominated for typical sparse in-vitro MEA graphs.
        %   Network: Density/GlobalEfficiency/Modularity are the most stable
        %     network-level descriptors; StdFR/RiseTime/FallTime are noisy.
            switch group
                case "ActivityFeatures"
                    names = ["FiringRate", "CV2InterSpikeInterval", ...
                             "RevisedLocalVariation", "FanoFactor"];

                case "WaveformFeatures"
                    names = ["T2Pdelay", "HalfWidth", "Asymmetry"];

                case "RegularityFeatures"
                    names = "RegularityFrequency";   % RegularityFit is NaN without a significant rhythm

                case "Catch22"
                    names = string.empty;   % excluded from core

                case "BombcellMetrics"
                    names = ["BC_fractionRPVs", "BC_waveformDuration", "BC_presenceRatio"];

                case "BurstParticipation"
                    names = unitFeatureNames("burst");

                case "GraphFeatures"
                    names = string.empty;   % excluded from core (too sparse / NaN-dominated)

                case "NetworkBurst"
                    names = ["MeanInterBurstInterval", "BurstDuration", "PeakFR", ...
                             "IntraFiringRate", "InterFiringRate"];

                case "NetworkRegularity"
                    names = "NetworkRegularityFrequency";

                case "NetworkCatch22"
                    names = string.empty;   % excluded from core

                case "NetworkGraph"
                    % Include stable network-level graph features only.
                    % These are column prefix-matched in selectFeatureCols because
                    % actual names are suffixed by algorithm (e.g. Density_CCG).
                    % Listing bare roots here; FeatureStore applies prefix matching.
                    names = ["Density", "GlobalEfficiency", "Modularity"];

                case "CellTypeGraph"
                    % Prefix-matched: Density_CCG_EE, Modularity_STTC_II, etc.
                    names = ["Density", "GlobalEfficiency", "Modularity", "CrossTypeDensity"];

                case "CellTypeBalance"
                    names = ["ExcitatoryFraction", "MeanFiringRate_E", ...
                             "MeanFiringRate_I", "FiringRateRatio_EI"];

                case "CellTypeActivity"
                    names = ["MeanCV2_E", "MeanCV2_I", ...
                             "MeanFanoFactor_E", "MeanFanoFactor_I"];

                case "CellTypeBurst"
                    % what the labels add beyond firing rate (computeCellTypeBurstEI)
                    names = ["BurstEI_IEContrastIn_Excess", "BurstEI_IFoldIn_Excess", ...
                             "BurstEI_EFoldIn_Excess", "BurstEI_IAfterPeak_Excess", ...
                             "BurstUnitLatency_IminusE_Excess", ...
                             "BurstUnitParticipation_IminusE_Excess"];

                case "CellTypeCorrelation"
                    names = ["MeanSTTC_EE", "MeanSTTC_II", "MeanSTTC_EI", ...
                             "STTCSynchronyIndex"];

                case "SpatialFeatures"
                    % Core spatial: extent + gradient + FC decay.
                    % E/I-specific (label-dependent), clustering (often NaN),
                    % per-unit, and burst propagation features excluded.
                    names = ["ConvexHullArea", "ChipCoverage", "MeanPairwiseDistance", ...
                             "SpatialFRMoransI", "CenterPeripheryFR_ratio", ...
                             "SpatialDecayTau", "DistanceFCCorrelation"];

                otherwise
                    names = string.empty;
            end
        end

        function names = fullSet(group)
        % FULLSET  All feature names for the given group (no filtering).
        %   Single-unit groups come from unitFeatureNames.
            switch group
                case "ActivityFeatures"
                    names = unitFeatureNames("act");

                case "WaveformFeatures"
                    % Actual column names written by inferWaveformFeatures.m
                    names = ["AUC_peak_1", "AUC_trough", "AUC_peak_2", ...
                             "Rise", "Decay", "HalfWidth", "Asymmetry", ...
                             "T2Pdelay", "T2Pratio", "Amplitude"];   % Amplitude: uV, NaN without a scale

                case "RegularityFeatures"
                    names = unitFeatureNames("reg");

                case "Catch22"
                    try
                        names = "SC_" + string(GetAllFeatureNames(false));   % catch22 (no DN_Mean / DN_Spread_Std)
                    catch
                        names = string.empty;
                    end

                case "BombcellMetrics"
                    names = unitFeatureNames("bc");

                case "BurstParticipation"
                    names = unitFeatureNames("burst");   % fraction of spikes in network bursts

                case "GraphFeatures"
                    % Per-unit graph columns: suffixed by algorithm (e.g. _CCG, _STTC)
                    % Listed as roots; actual names contain these as prefixes.
                    names = ["ClusteringCoefficient", "LocalEfficiency", "Strength_", ...
                             "EigenCentrality", "Betweenness"];

                case "NetworkBurst"
                    names = ["MeanInterBurstInterval", "BurstDuration", ...
                             "RiseTime", "FallTime", "PeakFR", "StdFR", ...
                             "StdBurstDuration", "StdInterBurstInterval", ...
                             "IntraFiringRate", "InterFiringRate", ...
                             "BurstWidth", "BurstRise", "BurstFall", "BurstFWHM"];

                case "NetworkRegularity"
                    names = "Network" + unitFeatureNames("reg");  % NetworkRegularityFrequency, ...

                case "NetworkCatch22"
                    names = string(GetAllFeatureNames(false));   % plain catch22 names (SC_ = single unit)

                case "NetworkGraph"
                    % Network-level graph columns: suffixed by algorithm
                    names = ["Density", "Assortativity", "RichClub", "MeanStrength_", ...
                             "GlobalEfficiency", "Modularity", "SmallWorldness"];

                case "CellTypeGraph"
                    % All graph roots (prefix-matched with _EE/_II suffixes)
                    % plus cross-type and per-unit connection counts
                    names = ["Density", "Assortativity", "RichClub", ...
                             "GlobalEfficiency", "Modularity", "SmallWorldness", ...
                             "CrossTypeDensity", "WithinTypeConnections", "CrossTypeConnections"];

                case "CellTypeBalance"
                    names = ["ExcitatoryFraction", "MeanFiringRate_E", ...
                             "MeanFiringRate_I", "FiringRateRatio_EI"];

                case "CellTypeActivity"
                    names = ["MeanCV2_E", "MeanCV2_I", ...
                             "MeanFanoFactor_E", "MeanFanoFactor_I", ...
                             "MeanLvR_E", "MeanLvR_I", ...
                             "MeanISIBimodality_E", "MeanISIBimodality_I"];

                case "CellTypeBurst"
                    [net_names, unit_names] = cellTypeBurstFeatureNames();
                    names = [net_names, unit_names];

                case "CellTypeCorrelation"
                    names = ["MeanSTTC_EE", "MeanSTTC_II", "MeanSTTC_EI", ...
                             "STTCSynchronyIndex", ...
                             "MeanSTTC_WithinType", "MeanSTTC_CrossType"];

                case "SpatialFeatures"
                    names = ["ConvexHullArea", "ChipCoverage", ...
                             "MeanPairwiseDistance", "CentroidSpread", ...
                             "SpatialMixingIndex", "MeanFractionExcNN_E", "MeanFractionExcNN_I", ...
                             "DistFromCentroid", "FractionExcNN", "FractionInhNN", ...
                             "SpatialFRMoransI", "CenterPeripheryFR_ratio", ...
                             "RipleysL_max", "ClusterScale", ...
                             "SpatialDecayTau", "DistanceFCCorrelation", "MeanFC_NearNeighbors", ...
                             "SpatialEIVariability", "SpatialEIMoransI", ...
                             "BurstOriginDispersion", "MeanBurstPropagationSpeed", "BurstOriginExcFraction"];

                otherwise
                    names = string.empty;
            end
        end

    end
end

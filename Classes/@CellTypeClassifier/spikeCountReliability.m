function R = spikeCountReliability(ctc, opts)
% SPIKECOUNTRELIABILITY  How many spikes a unit needs for a stable label (spike-train thinning).
%
%   R = ctc.spikeCountReliability()
%   R = ctc.spikeCountReliability(Counts=[50 100 200 500 1000], NUnits=500)
%
% Units with many spikes are randomly thinned to each count; their ACG is recomputed from the
% thinned train (waveform unchanged) and passed through the same smoothing, harmonisation,
% normalisation and composition model as in estimateComposition. The label is compared with
% the label of the same unit computed the same way from all its spikes. Thinning keeps the
% expected ACG shape and only adds the sampling noise of fewer spikes, so the agreement
% isolates the effect of the spike count from that of the (biologically different) low-rate
% units. Run estimateComposition first.
%
% OPTIONS:
%   Counts          - spike counts to test (default [25 50 100 200 400 800 1600 3200])
%   NUnits          - source units per count, drawn at random (default 1000)
%   MinSourceFactor - source units have >= MinSourceFactor x count spikes (default 4)
%   Units           - restrict the source units (indices into NormalizedFeatures.unique_ud),
%                     e.g. units with independent ground truth
%   TargetAgreement - MinSpikes = smallest count from which on the label agreement stays
%                     >= this (default 0.9)
%
% OUTPUT / SETS ctc.Composition.Reliability:
%   R.Table     - Count, NUnits, Agreement (label as with all spikes), AgreementConfident (among
%                 units with P outside [0.2 0.8] with all spikes), FlipToConfident (share of the
%                 thinned units with a confident (P outside [0.2 0.8]) WRONG label),
%                 MeanAbsDelta (|P - P_all|), Spearman (P vs P_all)
%   R.Units     - per thinned unit: UniqueUnit (index into NormalizedFeatures.unique_ud), Count,
%                 PAll, PThin (P(inhibitory) at the training prior, all vs thinned spikes)
%   R.MinSpikes - recommended minimum (NaN if no tested count reaches TargetAgreement);
%                 (apply with confidentLabels(MinSpikes=R.MinSpikes))
arguments
    ctc CellTypeClassifier
    opts.Counts double = [25 50 100 200 400 800 1600 3200]
    opts.NUnits (1,1) double = 1000
    opts.MinSourceFactor (1,1) double = 4
    opts.TargetAgreement (1,1) double = 0.9
    opts.Units double = []
end
assert(isfield(ctc.Composition, 'Model'), 'CellTypeClassifier:noComposition', ...
    'Run estimateComposition() before spikeCountReliability().');
nf = ctc.NormalizedFeatures;
ud = nf.unique_ud;
n_spk = arrayfun(@(u) numel(u.SpikeTimes), ud);
rs = RandStream('twister', 'Seed', ctc.Parameters.RNGSeed);
counts = sort(opts.Counts(:)');
nc = numel(counts);
[n_src, agree, agree_conf, flip_conf, mad, rho] = deal(nan(nc, 1));
units = cell(nc, 1);
for k = 1:nc
    src = find(n_spk >= opts.MinSourceFactor * counts(k));
    if ~isempty(opts.Units), src = intersect(src, opts.Units(:)'); end
    if numel(src) < 20, continue; end
    if numel(src) > opts.NUnits, src = src(sort(randperm(rs, numel(src), opts.NUnits))); end
    p_all  = scoreUnits(ctc, ud(src), src);
    thin = ud(src);
    for i = 1:numel(thin)
        thin(i).SpikeTimes = CellTypeClassifier.thinSpikeTrain(thin(i).SpikeTimes(:), counts(k), rs);
    end
    p_thin = scoreUnits(ctc, thin, src);
    same = (p_thin > 0.5) == (p_all > 0.5);
    conf = p_all >= 0.8 | p_all <= 0.2;
    n_src(k) = numel(src);
    agree(k) = mean(same);
    agree_conf(k) = mean(same(conf));
    flip_conf(k) = mean(~same & (p_thin >= 0.8 | p_thin <= 0.2));
    mad(k) = mean(abs(p_thin - p_all));
    rho(k) = corr(p_thin, p_all, 'Type', 'Spearman');
    units{k} = table(src(:), repmat(counts(k), numel(src), 1), p_all, p_thin, 'VariableNames', ...
        {'UniqueUnit', 'Count', 'PAll', 'PThin'});
end
T = table(counts(:), n_src, agree, agree_conf, flip_conf, mad, rho, 'VariableNames', ...
    {'Count', 'NUnits', 'Agreement', 'AgreementConfident', 'FlipToConfident', 'MeanAbsDelta', 'Spearman'});
% smallest tested count from which on every larger tested count also reaches the target
ok = T.Agreement >= opts.TargetAgreement | isnan(T.Agreement);
first = find(~ok, 1, 'last');
if isempty(first), first = 0; end
min_spikes = NaN;
if first < nc && ~isnan(T.Agreement(first + 1)), min_spikes = T.Count(first + 1); end
R = struct('Table', T, 'MinSpikes', min_spikes, 'TargetAgreement', opts.TargetAgreement, ...
    'Units', vertcat(units{:}));
ctc.Composition.Reliability = R;
disp(T);
fprintf('spikeCountReliability: label agreement >= %.0f%% from %g spikes on (%.0f%% of units have fewer)\n', ...
    100 * opts.TargetAgreement, min_spikes, 100 * mean(n_spk < min_spikes));
end

function p = scoreUnits(ctc, ud, unique_idx)
% P(inhibitory) of the composition model for UnitData ud (unique units unique_idx, for their
% per-group z-score statistics), computed as in buildNormalizedFeatures + estimateComposition.
nf = ctc.NormalizedFeatures; np = ctc.NormalizationParams; M = ctc.Composition.Model;
ph = ctc.Parameters.Harmonization;
N = numel(ud);
wf_list = {ud.ReferenceWaveform};
L = cellfun(@numel, wf_list);
wf = nan(max(L), N);
for i = 1:N, wf(1:L(i), i) = double(wf_list{i}(:)); end
sr = [ud.SamplingRate];
acg = zeros(round(ph.ACGLag / ph.ACGBinSize) * 2 + 1, N);
for i = 1:N
    acg(:, i) = CellTypeClassifier.computeACG(ud(i).SpikeTimes(:), ph.ACGBinSize, ph.ACGLag);
end
acg = smoothACG(acg, nf.acg_sigma_ms / 1000 / ph.ACGBinSize);
[~, X] = evalc('buildFeatureMatrix(ctc, wf, acg, sr)');
X(isnan(X)) = 0;
if isfield(nf, 'group_stats') && ~isempty(nf.group_stats)
    X = zscoreByGroup(X, nf.group_values(unique_idx), Stats=nf.group_stats, ConstantAs=NaN);
end
X = normalize(X, 'center', np.mu_global, 'scale', np.sigma_global);
X(:, np.nan_cols) = [];
X = X ./ np.scale;
if isfield(np, 'feature_selection_mask') && ~all(np.feature_selection_mask)
    X = X(:, np.feature_selection_mask);
end
X(~isfinite(X)) = 0;
[~, s] = predict(M.Classifier, (X - M.Mu) * M.Coeff);
p = s(:, 2);
end

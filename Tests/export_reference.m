function export_reference(outDir)
% EXPORT_REFERENCE  Write golden reference files of deterministic DeePhys functions.
%
%   export_reference(outDir)
%
% Runs deterministic functions on small synthetic inputs (fixed seeds) and stores inputs and
% outputs, one HDF5 file per case, to check another implementation (e.g. a Python port)
% numerically against this one:
%
%   outDir/<function>/<case>.h5        datasets /inputs/<name>, /outputs/<name>
%   outDir/<function>/<case>_<name>.parquet   tables
%   outDir/SOURCE.json                 repository commit, MATLAB / toolbox versions, case list
%   outDir/README.md                   file format description
%
% All unit / template / channel / community / seed indices are written 0-based (see README.md,
% and the root attribute "index_valued" of each file). Arrays are written as MATLAB holds them
% (column-major); readers with row-major arrays must reverse the axes (root attribute
% "layout"). Stochastic functions (Louvain communities, surrogate-based STTC significance,
% degree-preserving rewiring) are not exported.
%
% EXAMPLE:
%   export_reference(fullfile(tempdir, 'golden'))

    arguments
        outDir (1,1) string
    end

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root, 'Functions')));
    addpath(genpath(fullfile(root, 'Classes')));
    addpath(genpath(fullfile(root, 'Toolboxes')));
    addpath(fullfile(root, 'Tests', 'fixtures'));

    if ~isfolder(outDir), mkdir(outDir); end
    previous = cd(char(outDir));   % absolute path without Java
    outDir = string(pwd);
    cd(previous);

    % Mock Kilosort-style sorting (fixed seed inside make_mock_sorting), removed afterwards
    mock_parent = tempname;
    mkdir(mock_parent);
    cleanup = onCleanup(@() rmdir(mock_parent, 's'));
    n_mock = 20; fs_mock = 20000; dur_mock = 60;
    [~, mock_dir] = evalc('make_mock_sorting(mock_parent, n_mock, 150, fs_mock, dur_mock)');
    mock_dir = char(mock_dir);

    ctx = struct('outDir', outDir, 'cases', {{}}, 'log', containers.Map());
    ctx.log('cases') = {};

    [trains, dur, extra] = buildTrains(mock_dir, n_mock, fs_mock, dur_mock);

    ctx = exportCCG(ctx, trains);
    ctx = exportSTTC(ctx, trains, dur);
    ctx = exportActivity(ctx, trains, dur);
    ctx = exportRegularity(ctx, extra);
    ctx = exportWaveform(ctx, mock_dir);
    ctx = exportBursts(ctx);
    ctx = exportBCT(ctx);
    ctx = exportGraph(ctx);
    ctx = exportNormalization(ctx);
    ctx = exportACG(ctx, trains, fs_mock);

    writeSource(ctx, root);
    writeReadme(ctx.outDir);
    fprintf('export_reference: %d cases written to %s\n', numel(ctx.log('cases')), ctx.outDir);
end

% ======================================================================= inputs ===

function [trains, dur, extra] = buildTrains(mock_dir, n_mock, fs, dur)
% Spike trains: units of the mock sorting plus four synthetic ones (regular, bursty,
% rhythmic, two spikes). Unit order = unit id (1-based inside MATLAB).
    st  = double(readNPY(fullfile(mock_dir, 'spike_times.npy'))) / fs;
    tid = double(readNPY(fullfile(mock_dir, 'spike_templates.npy')));
    trains = cell(1, n_mock + 4);
    for u = 1:n_mock
        trains{u} = sort(st(tid == u - 1));
    end
    rng(11, 'twister');
    trains{n_mock + 1} = sort((0.5:0.5:dur - 0.5)' + 0.02 * randn(2 * dur - 1, 1));
    rng(12, 'twister');
    centres = (3:3:dur - 3)' + 0.3 * rand(dur / 3 - 1, 1);
    b = centres + 0.03 * rand(numel(centres), 8);
    trains{n_mock + 2} = sort(b(:));
    rng(13, 'twister');
    cand = sort(rand(round(dur * 12), 1) * dur);
    keep = rand(size(cand)) < (1 + 0.9 * sin(2 * pi * 2 * cand)) / 1.9;
    trains{n_mock + 3} = cand(keep);
    trains{n_mock + 4} = [10.2; 33.7];
    extra.rhythmic = trains{n_mock + 3};
    extra.duration = dur;
end

function [st, su] = mergeTrains(trains)
    st = vertcat(trains{:});
    su = cell2mat(arrayfun(@(u) u * ones(numel(trains{u}), 1), 1:numel(trains), 'UniformOutput', false)');
    [st, o] = sort(st);
    su = su(o);
end

function [st, su] = correlatedTrio(seed, T, n)
% Unit 2 follows unit 1 at 2 ms (jitter), unit 3 independent.
    rng(seed, 'twister');
    st1 = sort(rand(n, 1) * T);
    st2 = st1 + 0.002 + randn(n, 1) * 0.0005;
    st2 = st2(st2 > 0 & st2 < T);
    st3 = sort(rand(n, 1) * T);
    st = [st1; st2; st3];
    su = [ones(n, 1); 2 * ones(numel(st2), 1); 3 * ones(n, 1)];
    [st, o] = sort(st);
    su = su(o);
end

function [st, su] = suppressedPair(seed, T, n)
% Unit 2 is silenced 1-5 ms after each unit-1 spike (one-sided inhibition 1 -> 2).
    rng(seed, 'twister');
    st1 = sort(rand(n, 1) * T);
    st2 = sort(rand(n, 1) * T);
    d = st2 - st1';
    st2 = st2(~any(d > 0.001 & d <= 0.005, 2));
    st = [st1; st2];
    su = [ones(numel(st1), 1); 2 * ones(numel(st2), 1)];
    [st, o] = sort(st);
    su = su(o);
end

% ======================================================================= cases ====

function ctx = exportCCG(ctx, trains)
    [st, su] = correlatedTrio(21, 60, 1500);
    fs = 30000;
    in = struct('spike_times', st, 'spike_units', su);
    [ccg, t] = CCG(st, su, 'binSize', 0.001, 'duration', 0.1, 'Fs', 1/fs);
    ctx = putCase(ctx, 'CCG', 'trio_bin1ms_dur100ms_counts', ...
        struct('binSize', 0.001, 'duration', 0.1, 'Fs', 1/fs, 'norm', 'counts'), ...
        in, struct('ccg', ccg, 't', t), "inputs.spike_units");
    [ccg2, t2] = CCG(st, su, 'binSize', 0.002, 'duration', 0.2, 'Fs', 1/fs, 'norm', 'rate');
    ctx = putCase(ctx, 'CCG', 'trio_bin2ms_dur200ms_rate', ...
        struct('binSize', 0.002, 'duration', 0.2, 'Fs', 1/fs, 'norm', 'rate'), ...
        in, struct('ccg', ccg2, 't', t2), "inputs.spike_units");
    [mst, msu] = mergeTrains(trains(1:20));
    [ccg3, t3] = CCG(mst, msu, 'binSize', 0.005, 'duration', 1, 'Fs', 1/20000);
    ctx = putCase(ctx, 'CCG', 'mock_units_bin5ms_dur1s_counts', ...
        struct('binSize', 0.005, 'duration', 1, 'Fs', 1/20000, 'norm', 'counts'), ...
        struct('spike_times', mst, 'spike_units', msu), struct('ccg', ccg3, 't', t3), "inputs.spike_units");

    % bz_cch_conv on correlograms of the trio (columns: 1->2, 2->1, 1->3)
    cch = double([ccg(:, 1, 2), ccg(:, 2, 1), ccg(:, 1, 3)]);
    cases = {
        'gaussian_w10',       10, [],         [];
        'gaussian_w5_hf0p5',   5, 'gaussian', 0.5;
        'rect_w5',             5, 'rect',     [];
        'triang_w5',           5, 'triang',   []};
    for k = 1:size(cases, 1)
        w = cases{k, 2}; wt = cases{k, 3}; hf = cases{k, 4};
        if isempty(wt)
            [p, pred, q] = bz_cch_conv(cch, w);
        else
            [p, pred, q] = bz_cch_conv(cch, w, wt, hf);
        end
        opts = struct('W', w, 'WINTYPE', string(wt), 'HF', hf);
        ctx = putCase(ctx, 'bz_cch_conv', cases{k, 1}, opts, struct('CCH', cch), ...
            struct('pvals', p, 'pred', pred, 'qvals', q), strings(0));
    end
    [p, pred, q] = bz_cch_conv(cch(:, 1), 10);
    ctx = putCase(ctx, 'bz_cch_conv', 'single_column_gaussian_w10', struct('W', 10), ...
        struct('CCH', cch(:, 1)), struct('pvals', p, 'pred', pred, 'qvals', q), strings(0));

    % computeConnectivityCCG (CCG + convolution predictor + monosynaptic call)
    pr = struct('BinSize', 0.001, 'Duration', 0.1, 'Conv_w', 10, 'Alpha', 0.001);
    idxc = ["inputs.spike_units", "outputs.ExcitatoryConnection", "outputs.InhibitoryConnection"];
    for crit = ["strict", "cellexplorer"]
        p2 = pr; p2.InhibitoryCriterion = crit;
        [~, conn] = evalc('computeConnectivityCCG(st, su, 3, fs, p2)');
        o = struct('ExcitatoryConnection', conn.ExcitatoryConnection, ...
            'InhibitoryConnection', conn.InhibitoryConnection, 'bd', conn.bd, 'CCGs', conn.CCGs);
        opts = p2; opts.n_units = 3; opts.sampling_rate = fs;
        ctx = putCase(ctx, 'computeConnectivityCCG', "trio_" + crit, opts, ...
            struct('spike_times', st, 'spike_units', su), o, idxc);
    end
    [st2, su2] = suppressedPair(22, 300, 3000);
    p2 = pr; p2.InhibitoryCriterion = "strict";
    [~, conn] = evalc('computeConnectivityCCG(st2, su2, 2, 10000, p2)');
    o = struct('ExcitatoryConnection', conn.ExcitatoryConnection, ...
        'InhibitoryConnection', conn.InhibitoryConnection, 'bd', conn.bd, 'CCGs', conn.CCGs);
    opts = p2; opts.n_units = 2; opts.sampling_rate = 10000;
    ctx = putCase(ctx, 'computeConnectivityCCG', 'suppressed_pair_strict', opts, ...
        struct('spike_times', st2, 'spike_units', su2), o, idxc);
end

function ctx = exportSTTC(ctx, trains, dur)
    [st, su] = correlatedTrio(31, 60, 400);
    pair = @(a, b) struct('spike_times_1', st(su == a), 'spike_times_2', st(su == b), ...
        'N1', nnz(su == a), 'N2', nnz(su == b), 'Time', [0 60]);
    specs = {'corr_dt50ms', 1, 2, 0.05; 'corr_dt10ms', 1, 2, 0.01; 'indep_dt50ms', 1, 3, 0.05;
             'corr_dt2ms', 1, 2, 0.002};
    for k = 1:size(specs, 1)
        in = pair(specs{k, 2}, specs{k, 3});
        in.dt = specs{k, 4};
        v = sttc(in.N1, in.N2, in.dt, in.Time, in.spike_times_1, in.spike_times_2);
        ctx = putCase(ctx, 'sttc', specs{k, 1}, struct('dt', in.dt), in, struct('tileCoef', v), strings(0));
    end
    % degenerate pairs (NaN / single spike)
    in = struct('spike_times_1', 5, 'spike_times_2', [4.99; 20], 'N1', 1, 'N2', 2, 'Time', [0 60], 'dt', 0.05);
    ctx = putCase(ctx, 'sttc', 'single_spike_vs_two', struct('dt', 0.05), in, ...
        struct('tileCoef', sttc(1, 2, 0.05, [0 60], 5, [4.99; 20])), strings(0));
    in = struct('spike_times_1', zeros(0, 1), 'spike_times_2', [1; 2; 3], 'N1', 0, 'N2', 3, 'Time', [0 60], 'dt', 0.05);
    ctx = putCase(ctx, 'sttc', 'empty_train_nan', struct('dt', 0.05), in, ...
        struct('tileCoef', sttc(0, 3, 0.05, [0 60], zeros(0, 1), [1; 2; 3])), strings(0));
    % all pairs of the mock-sorting units plus synthetic units (last unit has no spikes)
    tr = [trains, {zeros(0, 1)}];
    n = numel(tr);
    for dt = [0.05 0.01]
        M = zeros(n);
        for i = 1:n - 1
            for j = i + 1:n
                M(i, j) = sttc(numel(tr{i}), numel(tr{j}), dt, [0 dur], tr{i}, tr{j});
            end
        end
        M = M + M';
        [mst, msu] = mergeTrains(tr);
        in = struct('spike_times', mst, 'spike_units', msu, 'duration', dur, 'n_units', n, 'dt', dt);
        ctx = putCase(ctx, 'sttc', sprintf('all_pairs_dt%dms', round(dt * 1000)), ...
            struct('dt', dt, 'note', 'symmetric matrix of sttc over unit pairs; diagonal 0'), ...
            in, struct('sttc_matrix', M), "inputs.spike_units");
    end
end

function ctx = exportActivity(ctx, trains, dur)
    [st, su] = mergeTrains(trains);
    n = numel(trains);
    runs = {'default', struct('FanoBinWidth', 0.1, 'RefractoryPeriod', 0.005);
            'fano1s_refr2ms_regularity', struct('FanoBinWidth', 1, 'RefractoryPeriod', 0.002, ...
                'Regularity', struct('Binning', 0.05, 'N_peaks', 3))};
    for r = 1:size(runs, 1)
        p = runs{r, 2};
        A = table(); R = table();
        for u = 1:n
            rng(100 + u, 'twister');   % the decay fit (curve fitting toolbox) may start from a random point
            [a, rg] = computeActivityFeatures(trains{u}, dur, p);
            A = [A; a]; %#ok<AGROW>
            if ~isempty(rg), R = [R; rg]; end %#ok<AGROW>
        end
        o = struct('activity', A);
        if height(R) > 0, o.regularity = R; end
        ctx = putCase(ctx, 'computeActivityFeatures', runs{r, 1}, p, ...
            struct('spike_times', st, 'spike_units', su, 'duration', dur, 'n_units', n), o, "inputs.spike_units");
    end
    % optional catch22 block (needs the catch22 mex files)
    try
        p = struct('FanoBinWidth', 0.1, 'RefractoryPeriod', 0.005, 'Catch22', struct('BinSize', 0.05));
        C = table();
        for u = 1:n
            [~, ~, c22] = computeActivityFeatures(trains{u}, dur, p);
            C = [C; c22]; %#ok<AGROW>
        end
        ctx = putCase(ctx, 'computeActivityFeatures', 'catch22_bin50ms', p, ...
            struct('spike_times', st, 'spike_units', su, 'duration', dur, 'n_units', n), ...
            struct('catch22', C), "inputs.spike_units");
    catch ME
        warning('export_reference:catch22', 'catch22 case skipped: %s', ME.message);
    end
end

function ctx = exportRegularity(ctx, extra)
    dur = extra.duration;
    binned = @(st, bw) histcounts(st, 'BinLimits', [0 dur], 'BinWidth', bw);
    norm1 = @(x) x / max(x);
    s1 = norm1(binned(extra.rhythmic, 0.05));
    cases = {'rhythmic_train_bin50ms_peaks3', s1, 0.05, 3};
    cases(end+1, :) = {'rhythmic_train_bin50ms_peaks5', s1, 0.05, 5};
    rng(41, 'twister');
    poisson_train = sort(rand(600, 1) * dur);
    cases(end+1, :) = {'poisson_train_bin100ms_peaks3', norm1(binned(poisson_train, 0.1)), 0.1, 3};
    rng(42, 'twister');
    t = (0:0.01:10 - 0.01);
    sig = sin(2 * pi * 3 * t) + 0.5 * sin(2 * pi * 6 * t) + 0.25 * sin(2 * pi * 9 * t) + 0.3 * randn(size(t));
    cases(end+1, :) = {'sine_with_harmonics_bin10ms_peaks3', sig, 0.01, 3};
    for k = 1:size(cases, 1)
        rng(200 + k, 'twister');   % the decay fit (curve fitting toolbox) may start from a random point
        [f, m, c] = computeRegularity(cases{k, 2}, cases{k, 3}, cases{k, 4});
        ctx = putCase(ctx, 'computeRegularity', cases{k, 1}, ...
            struct('binning', cases{k, 3}, 'n_peaks', cases{k, 4}, ...
                'note', 'fit_coeff comes from a nonlinear fit that may start at a random point; compare it with a loose tolerance'), ...
            struct('binned_signal', cases{k, 2}), ...
            struct('frequency', f, 'magnitude', m, 'fit_coeff', c), strings(0));
    end
end

function ctx = exportWaveform(ctx, mock_dir)
    tm = double(readNPY(fullfile(mock_dir, 'templates.npy')));   % units x samples x channels
    nu = size(tm, 1);
    wf = zeros(size(tm, 2), nu);
    pc = zeros(1, nu);
    for u = 1:nu
        [~, pc(u)] = min(min(squeeze(tm(u, :, :)), [], 1));
        w = squeeze(tm(u, :, pc(u)))';
        wf(:, u) = w / abs(min(w));
    end
    ctx = putCase(ctx, 'computeWaveformFeatures', 'mock_templates_fs20000', struct('sampling_rate', 20000), ...
        struct('norm_wf_matrix', wf, 'templates', tm, 'peak_channel', pc), ...
        struct('features', wfFeatureMatrix(computeWaveformFeatures(wf, 20000)), 'feature_names', wfNames()), ...
        "inputs.peak_channel");
    % parametric biphasic waveforms of different width at 30 kHz
    fs = 30000;
    t = (0:89)' / fs * 1000;   % ms
    w2 = zeros(numel(t), 6);
    for k = 1:6
        wd = 0.08 + 0.04 * k;
        w = -exp(-(t - 1.0).^2 / (2 * wd^2)) + (0.15 + 0.05 * k) * exp(-(t - 1.0 - 4 * wd).^2 / (2 * (3 * wd)^2)) ...
            + 0.1 * exp(-(t - 1.0 + 3 * wd).^2 / (2 * wd^2));
        w2(:, k) = w / abs(min(w));
    end
    ctx = putCase(ctx, 'computeWaveformFeatures', 'parametric_fs30000', struct('sampling_rate', fs), ...
        struct('norm_wf_matrix', w2), ...
        struct('features', wfFeatureMatrix(computeWaveformFeatures(w2, fs)), 'feature_names', wfNames()), strings(0));
    t = linspace(0, 3, 60)';
    w3 = 0.3 * exp(-((t - 1.0).^2) / 0.02) - exp(-((t - 1.3).^2) / 0.03) + 0.4 * exp(-((t - 1.8).^2) / 0.05);
    w3 = w3 / max(abs(w3));
    ctx = putCase(ctx, 'computeWaveformFeatures', 'single_unit_fs10000', struct('sampling_rate', 10000), ...
        struct('norm_wf_matrix', w3), ...
        struct('features', wfFeatureMatrix(computeWaveformFeatures(w3, 10000)), 'feature_names', wfNames()), strings(0));
end

function M = wfFeatureMatrix(feats)
    names = wfNames();
    M = nan(numel(feats), numel(names));
    for u = 1:numel(feats)
        for k = 1:numel(names)
            M(u, k) = feats{u}.(names(k));
        end
    end
end

function n = wfNames()
    n = ["AUC_peak_1", "AUC_trough", "AUC_peak_2", "Rise", "Decay", "HalfWidth", "Asymmetry", "T2Pdelay", "T2Pratio"];
end

function ctx = exportBursts(ctx)
    % burstShape on simulated bursts with known spans
    [st, su, truth] = simulateNetworkBursts(NUnits=30, Duration=120, Seed=2);
    dur = 120;
    shapeRuns = {'default_level0p25', struct(), truth.T_start, truth.T_end;
                 'level0p5_wide_kernel', struct('Level', 0.5, 'Bin', 0.005, 'Sigma', 0.01, 'Pad', 0.1), truth.T_start, truth.T_end;
                 'narrow_spans_level0p2', struct('Level', 0.2), truth.T_start + 0.05, truth.T_end - 0.05};
    for k = 1:size(shapeRuns, 1)
        a = namedargs2cell(shapeRuns{k, 2});
        S = burstShape(st, shapeRuns{k, 3}, shapeRuns{k, 4}, dur, a{:});
        ctx = putCase(ctx, 'burstShape', shapeRuns{k, 1}, shapeRuns{k, 2}, ...
            struct('spike_times', st, 'T_start', shapeRuns{k, 3}, 'T_end', shapeRuns{k, 4}, 'duration', dur), ...
            struct('shape', S), strings(0));
    end

    % detectNetworkBursts (Seed fixed; deterministic)
    sims = {'network', {'NUnits', 30, 'Duration', 120, 'Seed', 1}, {};
            'network_defaults_seed3', {'NUnits', 30, 'Duration', 120, 'Seed', 3}, {};
            'fixed_N_threshold', {'NUnits', 30, 'Duration', 120, 'Seed', 1}, {'N', 10, 'Threshold', 0.1, 'TonicMethod', 'none'};
            'coupling_tonic_10_surrogates', {'NUnits', 30, 'Duration', 120, 'Seed', 4}, {'TonicMethod', 'coupling', 'NSurrogates', 10};
            'nested_subbursts', {'NUnits', 30, 'Duration', 120, 'SubBursts', 3, 'BurstDuration', 0.6, 'SubDutyCycle', 0.7, 'Seed', 1}, {};
            'no_network_bursts_control', {'NUnits', 30, 'Duration', 120, 'Mode', 'none', 'Seed', 5}, {}};
    for k = 1:size(sims, 1)
        [st, su] = simulateNetworkBursts(sims{k, 2}{:});
        dur = sims{k, 2}{find(strcmp(sims{k, 2}, 'Duration')) + 1};
        optargs = [sims{k, 3}, {'Seed', 1}];
        [b, info] = detectNetworkBursts(st, su, dur, optargs{:});
        o = struct('bursts', b, 'info', info);
        opts = struct('simulation', {sims{k, 2}}, 'detector', {optargs});
        ctx = putCase(ctx, 'detectNetworkBursts', sims{k, 1}, opts, ...
            struct('spike_times', st, 'spike_units', su, 'duration', dur), o, ...
            ["inputs.spike_units", "outputs.bursts.Superburst", "outputs.info.Units.Unit"]);
    end
end

function ctx = exportBCT(ctx)
    n = 24;
    rng(101, 'twister');
    Bd = double(rand(n) < 0.18); Bd(1:n+1:end) = 0;
    U = triu(rand(n) < 0.18, 1); Bu = double(U | U');
    rng(102, 'twister');
    Wd = Bd .* (0.1 + rand(n)); Wd = Wd / max(Wd(:));
    Wu0 = triu(Bu .* (0.1 + rand(n)), 1); Wu = Wu0 + Wu0'; Wu = Wu / max(Wu(:));
    mats = struct('random_n24', struct('Bd', Bd, 'Bu', Bu, 'Wd', Wd, 'Wu', Wu));
    iso = @(A) setIso(A);
    mats.isolated_node_n24 = struct('Bd', iso(Bd), 'Bu', iso(Bu), 'Wd', iso(Wd), 'Wu', iso(Wu));
    names = fieldnames(mats);
    for m = 1:numel(names)
        c = names{m}; M = mats.(c);
        [k1, n1, e1] = density_dir(M.Bd);
        ctx = putCase(ctx, 'density_dir', c, struct(), struct('A', M.Bd), struct('kden', k1, 'N', n1, 'K', e1), strings(0));
        [k1, n1, e1] = density_und(M.Bu);
        ctx = putCase(ctx, 'density_und', c, struct(), struct('A', M.Bu), struct('kden', k1, 'N', n1, 'K', e1), strings(0));
        ctx = putCase(ctx, 'clustering_coef_bd', c, struct(), struct('A', M.Bd), struct('C', clustering_coef_bd(M.Bd)), strings(0));
        ctx = putCase(ctx, 'clustering_coef_bu', c, struct(), struct('A', M.Bu), struct('C', clustering_coef_bu(M.Bu)), strings(0));
        ctx = putCase(ctx, 'clustering_coef_wd', c, struct(), struct('W', M.Wd), struct('C', clustering_coef_wd(M.Wd)), strings(0));
        ctx = putCase(ctx, 'clustering_coef_wu', c, struct(), struct('W', M.Wu), struct('C', clustering_coef_wu(M.Wu)), strings(0));
        ctx = putCase(ctx, 'assortativity_bin', c + "_directed", struct('flag', 1), struct('A', M.Bd), ...
            struct('r', assortativity_bin(M.Bd, 1)), strings(0));
        ctx = putCase(ctx, 'assortativity_bin', c + "_undirected", struct('flag', 0), struct('A', M.Bu), ...
            struct('r', assortativity_bin(M.Bu, 0)), strings(0));
        ctx = putCase(ctx, 'assortativity_wei', c + "_directed", struct('flag', 1), struct('W', M.Wd), ...
            struct('r', assortativity_wei(M.Wd, 1)), strings(0));
        ctx = putCase(ctx, 'assortativity_wei', c + "_undirected", struct('flag', 0), struct('W', M.Wu), ...
            struct('r', assortativity_wei(M.Wu, 0)), strings(0));
        [R, Nk, Ek] = rich_club_bd(M.Bd, 5);
        ctx = putCase(ctx, 'rich_club_bd', c, struct('k_level', 5), struct('A', M.Bd), struct('R', R, 'Nk', Nk, 'Ek', Ek), strings(0));
        [R, Nk, Ek] = rich_club_bu(M.Bu, 5);
        ctx = putCase(ctx, 'rich_club_bu', c, struct('k_level', 5), struct('A', M.Bu), struct('R', R, 'Nk', Nk, 'Ek', Ek), strings(0));
        ctx = putCase(ctx, 'efficiency_bin', c + "_global_directed", struct('local', 0), struct('A', M.Bd), struct('E', efficiency_bin(M.Bd)), strings(0));
        ctx = putCase(ctx, 'efficiency_bin', c + "_local_directed", struct('local', 1), struct('A', M.Bd), struct('E', efficiency_bin(M.Bd, 1)), strings(0));
        ctx = putCase(ctx, 'efficiency_bin', c + "_global_undirected", struct('local', 0), struct('A', M.Bu), struct('E', efficiency_bin(M.Bu)), strings(0));
        ctx = putCase(ctx, 'efficiency_bin', c + "_local_undirected", struct('local', 1), struct('A', M.Bu), struct('E', efficiency_bin(M.Bu, 1)), strings(0));
        ctx = putCase(ctx, 'efficiency_wei', c + "_global_directed", struct('local', 0), struct('W', M.Wd), struct('E', efficiency_wei(M.Wd)), strings(0));
        ctx = putCase(ctx, 'efficiency_wei', c + "_global_undirected", struct('local', 0), struct('W', M.Wu), struct('E', efficiency_wei(M.Wu)), strings(0));
        ctx = putCase(ctx, 'efficiency_wei', c + "_local2_undirected", struct('local', 2), struct('W', M.Wu), struct('E', efficiency_wei(M.Wu, 2)), strings(0));
        ctx = putCase(ctx, 'eigenvector_centrality_dir', c, struct(), struct('A', M.Bd), struct('v', eigenvector_centrality_dir(M.Bd)), strings(0));
        ctx = putCase(ctx, 'eigenvector_centrality_und', c, struct(), struct('A', M.Bu), struct('v', eigenvector_centrality_und(M.Bu)), strings(0));
        ctx = putCase(ctx, 'betweenness_bin', c + "_directed", struct(), struct('A', M.Bd), struct('BC', betweenness_bin(M.Bd)), strings(0));
        ctx = putCase(ctx, 'betweenness_bin', c + "_undirected", struct(), struct('A', M.Bu), struct('BC', betweenness_bin(M.Bu)), strings(0));
        Ld = weight_conversion(M.Wd, 'lengths'); Lu = weight_conversion(M.Wu, 'lengths');
        ctx = putCase(ctx, 'weight_conversion_lengths', c + "_directed", struct(), struct('W', M.Wd), struct('L', Ld), strings(0));
        ctx = putCase(ctx, 'betweenness_wei', c + "_directed", struct('note', 'input is weight_conversion(W, lengths)'), ...
            struct('L', Ld), struct('BC', betweenness_wei(Ld)), strings(0));
        ctx = putCase(ctx, 'betweenness_wei', c + "_undirected", struct('note', 'input is weight_conversion(W, lengths)'), ...
            struct('L', Lu), struct('BC', betweenness_wei(Lu)), strings(0));
        [is, os, st] = strengths_dir(M.Wd);
        ctx = putCase(ctx, 'strengths_dir', c, struct(), struct('W', M.Wd), struct('is', is, 'os', os, 'str', st), strings(0));
        ctx = putCase(ctx, 'strengths_und', c, struct(), struct('W', M.Wu), struct('str', strengths_und(M.Wu)), strings(0));
        D = distance_bin(M.Bd);
        ctx = putCase(ctx, 'distance_bin', c + "_directed", struct(), struct('A', M.Bd), struct('D', D), strings(0));
        [lam, eff, ecc, rad, dia] = charpath(D, 0, 0);
        ctx = putCase(ctx, 'charpath', c + "_finite_only", struct('diagonal_dist', 0, 'infinite_dist', 0), struct('D', D), ...
            struct('lambda', lam, 'efficiency', eff, 'ecc', ecc, 'radius', rad, 'diameter', dia), strings(0));
    end
end

function A = setIso(A)
    A(end, :) = 0; A(:, end) = 0;
end

function ctx = exportGraph(ctx)
    % two well separated clusters plus a small third group, 6 features
    rng(51, 'twister');
    X = [randn(30, 6); randn(30, 6) + 8; randn(15, 6) - 8];
    for k = [10 5]
        W = full(CellTypeClassifier.knnGraph(X, k));
        ctx = putCase(ctx, 'buildKNNGraph', sprintf('three_clusters_k%d', k), struct('k', k), ...
            struct('X', X), struct('W', W), strings(0));
    end
    Xd = [zeros(3, 2); 1 0; 0 1; 1 1];
    ctx = putCase(ctx, 'buildKNNGraph', 'duplicate_rows_k2', struct('k', 2), struct('X', Xd), ...
        struct('W', full(CellTypeClassifier.knnGraph(Xd, 2))), strings(0));
    rng(52, 'twister');
    Xs = randn(8, 3);
    ctx = putCase(ctx, 'buildKNNGraph', 'k_larger_than_n', struct('k', 100), struct('X', Xs), ...
        struct('W', full(CellTypeClassifier.knnGraph(Xs, 100))), strings(0));

    % label propagation on the first graph; the third cluster holds no seeds (unreached)
    G = full(CellTypeClassifier.knnGraph(X, 8));
    seedsets = {'balanced_seeds', [1:5, 31:35], [ones(1, 5), 2 * ones(1, 5)], 100, 1e-6;
                'unbalanced_seeds', [1:8, 31:33], [ones(1, 8), 2 * ones(1, 3)], 100, 1e-6;
                'not_converged_max_iter3', [1:5, 31:35], [ones(1, 5), 2 * ones(1, 5)], 3, 1e-9};
    for k = 1:size(seedsets, 1)
        [F, unreached, conv, n_iter] = CellTypeClassifier.propagateLabels(sparse(G), seedsets{k, 2}, seedsets{k, 3}, seedsets{k, 4}, seedsets{k, 5});
        ctx = putCase(ctx, 'propagateLabels', seedsets{k, 1}, ...
            struct('max_iter', seedsets{k, 4}, 'tol', seedsets{k, 5}), ...
            struct('G', G, 'seed_idx', seedsets{k, 2}, 'seed_y', seedsets{k, 3}), ...
            struct('F', F, 'unreached', unreached, 'converged', conv, 'n_iter', n_iter), "inputs.seed_idx");
    end

    % community enrichment: 6 communities, 3 cultures; community 2 holds responders, 5 few
    rng(61, 'twister');
    N = 150;
    M = randi(6, N, 1);
    culture = 10 * randi(3, N, 1);
    elig = rand(N, 1) < 0.9;
    p = 0.15 + 0.5 * (M == 2) - 0.12 * (M == 5) + 0.15 * (culture == 20);
    resp = rand(N, 1) < p;
    for mode = ["pooled", "stratified"]
        if mode == "pooled", s = []; else, s = culture; end
        [T, inh, exc] = CellTypeClassifier.communityEnrichment(M, elig, resp, 0.05, s);
        in = struct('M', M, 'eligible', elig, 'responsive', resp);
        if ~isempty(s), in.strata = s; end
        ctx = putCase(ctx, 'communityEnrichmentTest', mode, struct('fdr', 0.05, 'strata', mode), in, ...
            struct('T', T, 'inh_comm', inh, 'exc_comm', exc), ...
            ["inputs.M", "outputs.T.Community", "outputs.inh_comm", "outputs.exc_comm"]);
    end
end

function ctx = exportNormalization(ctx)
    step = {struct('type', 'combat', 'params', struct())};
    % 1. two batches, no covariate
    rng(71, 'twister');
    X = [randn(30, 5); randn(30, 5) * 1.5 + 4];
    batch = [ones(30, 1); 2 * ones(30, 1)];
    np = NormalizationPipeline(step);
    [Xo, np] = np.fit_transform(X, batch);
    Xt = X(1:10, :) + 0.1;
    ctx = putCase(ctx, 'ComBat', 'two_batches', struct('batch', 'numeric labels'), ...
        struct('X', X, 'batch', batch, 'X_new', Xt, 'batch_new', batch(1:10)), ...
        struct('X_out', Xo, 'X_new_out', np.transform(Xt, batch(1:10)), 'fit', combatFit(np)), strings(0));
    % 2. three string-labelled batches with a protected covariate
    rng(72, 'twister');
    n = [25 35 20];
    lab = repelem(["a"; "b"; "c"], n)';
    lab = lab(:);
    cov = [ones(15, 1); 2 * ones(10, 1); ones(10, 1); 2 * ones(25, 1); ones(12, 1); 2 * ones(8, 1)];
    X = randn(sum(n), 6);
    X(:, 1) = X(:, 1) + 2 * (cov == 1) - 2 * (cov == 2);
    X(lab == "b", :) = X(lab == "b", :) * 1.8 + 3;
    X(lab == "c", :) = X(lab == "c", :) - 2;
    np = NormalizationPipeline(step);
    [Xo, np] = np.fit_transform(X, lab, CovariateLabels=cov);
    ctx = putCase(ctx, 'ComBat', 'three_batches_covariate', struct('batch', 'string labels'), ...
        struct('X', X, 'batch', lab, 'covariate', cov), ...
        struct('X_out', Xo, 'fit', combatFit(np)), strings(0));
    % 3. missing values
    rng(73, 'twister');
    X = [randn(30, 5); randn(30, 5) + 3];
    X(randperm(300, 12)) = NaN;
    batch = [ones(30, 1); 2 * ones(30, 1)];
    np = NormalizationPipeline(step);
    [Xo, np] = np.fit_transform(X, batch);
    ctx = putCase(ctx, 'ComBat', 'missing_values', struct('batch', 'numeric labels'), ...
        struct('X', X, 'batch', batch), struct('X_out', Xo, 'fit', combatFit(np)), strings(0));
    % 4. transform data of an unseen batch (its location / scale come from its own rows)
    rng(74, 'twister');
    Xn = randn(12, 5) + 6;
    bn = 3 * ones(12, 1);
    ctx = putCase(ctx, 'ComBat', 'transform_unseen_batch', struct('batch', 'numeric labels'), ...
        struct('X', X, 'batch', batch, 'X_new', Xn, 'batch_new', bn), ...
        struct('X_out', Xo, 'X_new_out', np.transform(Xn, bn), 'fit', combatFit(np)), strings(0));

    % per-group z-scoring with reference rows
    rng(75, 'twister');
    X = [randn(30, 4) + 2; randn(30, 4) * 2 - 1; randn(6, 4) + 5];
    X(:, 4) = 3;                              % constant column
    g = repelem(["p1"; "p2"; "p3"], [30 30 6]);
    [Z, fb, stats] = zscoreByGroup(X, g);
    ctx = putCase(ctx, 'zscoreByGroup', 'all_rows', struct('min_ref', 20, 'ConstantAs', 0), ...
        struct('X', X, 'groups', g), struct('X', Z, 'fallback', fb, 'stats', stats), "outputs.fallback");
    is_ref = [true(15, 1); false(15, 1); true(10, 1); false(20, 1); true(6, 1)];
    [Z, fb, stats] = zscoreByGroup(X, g, is_ref, 12, ConstantAs=NaN);
    ctx = putCase(ctx, 'zscoreByGroup', 'reference_rows_constant_nan', struct('min_ref', 12, 'ConstantAs', NaN), ...
        struct('X', X, 'groups', g, 'is_ref', is_ref), struct('X', Z, 'fallback', fb, 'stats', stats), "outputs.fallback");
    Xn = randn(8, 4); gn = repelem(["p1"; "p2"], 4);
    Zn = zscoreByGroup(Xn, gn, Stats=stats, ConstantAs=NaN);
    ctx = putCase(ctx, 'zscoreByGroup', 'apply_stored_stats', struct('ConstantAs', NaN), ...
        struct('X', Xn, 'groups', gn, 'stats', stats), struct('X', Zn), strings(0));
end

function f = combatFit(np)
    c = np.FitParams.combat;
    f = struct('batches', string(c.batches), 'alpha', c.alpha, 'sigma', c.sigma, ...
        'gamma_star', c.gamma_star, 'delta_star', c.delta_star);
end

function ctx = exportACG(ctx, trains, fs)
    runs = {'unit1_bin0p5ms_lag100ms', 1, 0.0005, 0.1; 'unit3_bin1ms_lag50ms', 3, 0.001, 0.05;
            'rhythmic_unit_bin1ms_lag200ms', 23, 0.001, 0.2};
    for k = 1:size(runs, 1)
        st = trains{runs{k, 2}};
        acg = CellTypeClassifier.computeACG(st, runs{k, 3}, runs{k, 4}, fs);
        ctx = putCase(ctx, 'computeACG', runs{k, 1}, struct('bin_size', runs{k, 3}, 'lag', runs{k, 4}, 'sampling_rate', fs), ...
            struct('spike_times', st), struct('acg', acg), strings(0));
    end
end

% ======================================================================= writers ==

function ctx = putCase(ctx, fn, name, options, in, out, idx)
% Write one case file. idx: dot paths ("inputs.x", "outputs.T.col") of index-valued numeric
% data, held 1-based in MATLAB and written 0-based.
    fn = char(fn); name = char(name);
    d = fullfile(ctx.outDir, fn);
    if ~isfolder(d), mkdir(d); end
    f = fullfile(d, [name '.h5']);
    if isfile(f), delete(f); end
    idx = string(idx);
    for k = 1:numel(idx)
        parts = cellstr(split(idx(k), '.'));
        switch parts{1}
            case 'inputs',  in  = shiftPath(in,  parts(2:end));
            case 'outputs', out = shiftPath(out, parts(2:end));
            otherwise, error('export_reference:badIndexPath', 'Bad index path %s', idx(k));
        end
    end
    rec = struct('strings', containers.Map(), 'tables', containers.Map(), ...
        'empty', containers.Map(), 'lists', containers.Map());
    base = [fn '/' name];
    writeValue(f, '/inputs',  in,  rec, d, name);
    writeValue(f, '/outputs', out, rec, d, name);
    h5writeatt(f, '/', 'layout', 'matlab-column-major');
    h5writeatt(f, '/', 'function', fn);
    h5writeatt(f, '/', 'case', name);
    h5writeatt(f, '/', 'options', jsonencode(options));
    h5writeatt(f, '/', 'strings', jsonencode(rec.strings));
    h5writeatt(f, '/', 'tables', jsonencode(rec.tables));
    h5writeatt(f, '/', 'empty', jsonencode(rec.empty));
    h5writeatt(f, '/', 'lists', jsonencode(rec.lists));
    h5writeatt(f, '/', 'index_valued', jsonencode(cellstr(strrep(strrep(idx, '.', '/'), '//', '/'))));
    c = ctx.log('cases');
    c{end+1} = struct('function', fn, 'case', name, 'file', [base '.h5']);
    ctx.log('cases') = c;
end

function S = shiftPath(S, parts)
    if isempty(parts)
        if ~isnumeric(S), error('export_reference:badIndex', 'Index data must be numeric'); end
        S = S - 1;
    elseif (istable(S) && ~ismember(parts{1}, S.Properties.VariableNames)) || (isstruct(S) && ~isfield(S, parts{1}))
        return   % absent in this case (e.g. an empty result)
    elseif isstruct(S) || istable(S)
        S.(parts{1}) = shiftPath(S.(parts{1}), parts(2:end));
    else
        error('export_reference:badIndexPath', 'Cannot descend into %s', parts{1});
    end
end

function writeValue(f, p, v, rec, dirName, caseName)
    key = p(2:end);
    if istable(v)
        if width(v) == 0 || height(v) == 0
            rec.empty(key) = [height(v), width(v)];
            return
        end
        flat = strrep(regexprep(key, '^(inputs|outputs)/', ''), '/', '_');
        fname = [caseName '_' flat '.parquet'];
        pf = fullfile(dirName, fname);
        if isfile(pf), delete(pf); end
        parquetwrite(pf, v);
        rec.tables(key) = fname;
    elseif isstruct(v)
        if isscalar(v)
            fn = fieldnames(v);
            for k = 1:numel(fn)
                writeValue(f, [p '/' fn{k}], v.(fn{k}), rec, dirName, caseName);
            end
        else
            rec.lists(key) = numel(v);
            for k = 1:numel(v)
                writeValue(f, sprintf('%s/%d', p, k - 1), v(k), rec, dirName, caseName);
            end
        end
    elseif iscellstr(v) || ischar(v) || isstring(v) || iscategorical(v) %#ok<ISCLSTR>
        rec.strings(key) = jsonencode(cellstr(string(v)));
    elseif iscell(v)
        rec.lists(key) = numel(v);
        for k = 1:numel(v)
            writeValue(f, sprintf('%s/%d', p, k - 1), v{k}, rec, dirName, caseName);
        end
    elseif isnumeric(v) || islogical(v)
        if issparse(v), v = full(v); end
        if islogical(v), v = uint8(v); end
        if isempty(v)
            rec.empty(key) = size(v);
            return
        end
        h5create(f, p, size(v), 'Datatype', class(v));
        h5write(f, p, v);
    else
        error('export_reference:unsupportedType', 'Cannot write %s of class %s', key, class(v));
    end
end

function writeSource(ctx, root)
    [~, commit] = system(sprintf('git -C "%s" rev-parse HEAD', root));
    [~, desc]   = system(sprintf('git -C "%s" describe --tags --always --dirty', root));
    wanted = ["MATLAB", "Statistics and Machine Learning Toolbox", "Signal Processing Toolbox", ...
              "Curve Fitting Toolbox", "Econometrics Toolbox", "Parallel Computing Toolbox"];
    v = ver;
    tb = struct();
    for k = 1:numel(wanted)
        i = find(strcmp({v.Name}, wanted(k)), 1);
        if ~isempty(i), tb.(matlab.lang.makeValidName(wanted(k))) = v(i).Version; end
    end
    cases = ctx.log('cases');
    S = struct('deephys_commit', strtrim(commit), 'deephys_describe', strtrim(desc), ...
        'matlab_version', version, 'toolboxes', tb, ...
        'date', char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
        'n_cases', numel(cases), 'cases', {[cases{:}]});
    fid = fopen(fullfile(ctx.outDir, 'SOURCE.json'), 'w');
    fwrite(fid, jsonencode(S, 'PrettyPrint', true));
    fclose(fid);
end

function writeReadme(outDir)
    L = {
'# DeePhys golden reference files'
''
'Inputs and outputs of deterministic DeePhys (MATLAB) functions on small synthetic data,'
'for checking another implementation numerically. Written by `Tests/export_reference.m`.'
'`SOURCE.json` records the repository commit (`git describe`), MATLAB and toolbox versions,'
'the export date and the list of cases. Everything is derived from fixed seeds, so a second'
'export gives identical numbers (only the date in `SOURCE.json` differs).'
''
'## Layout'
''
'    <function>/<case>.h5                    one HDF5 file per case'
'    <function>/<case>_<name>.parquet        tables of that case'
''
'`<function>` is the MATLAB function (static wrappers: `buildKNNGraph` = `CellTypeClassifier.knnGraph`,'
'`propagateLabels`, `communityEnrichmentTest` = `CellTypeClassifier.communityEnrichment`,'
'`computeACG` = `CellTypeClassifier.computeACG`; `ComBat` = the `combat` step of `NormalizationPipeline`'
'(fit_transform / transform); `weight_conversion_lengths` = `weight_conversion(W, ''lengths'')`).'
''
'## HDF5 file'
''
'- Datasets `/inputs/<name>` and `/outputs/<name>`; structs become groups (`/outputs/info/Landscape/N`),'
'  cell arrays and struct arrays become groups with datasets `0`, `1`, ... (0-based).'
'- Arrays are stored as written by MATLAB `h5create` / `h5write` (column-major). h5py shows them'
'  with reversed axes: a MATLAB array of size [r c] is read as shape (c, r), a 3-D array [a b c] as'
'  (c, b, a). The root attribute `layout = "matlab-column-major"` marks this; a Python loader has to'
'  transpose every dataset (`a.T`, i.e. reverse all axes) to recover the MATLAB array. Vectors are'
'  therefore (1, n) or (n, 1) as in MATLAB; scalars are (1, 1).'
'- Logical values are `uint8` (0 / 1). NaN is stored as NaN. Numeric class is kept (double, single,'
'  uint16 for CCG counts, ...).'
'- Root attributes (all JSON strings except `layout`, `function`, `case`):'
'  - `options`      options / parameters the function was called with (NaN becomes null)'
'  - `strings`      map dataset path -> list of strings (strings and cellstr are not datasets)'
'  - `tables`       map dataset path -> parquet file name next to the h5 file'
'  - `empty`        map dataset path -> MATLAB size of an empty array / table that was not written'
'  - `lists`        map group path -> number of elements (cell arrays, struct arrays)'
'  - `index_valued` datasets / table columns that hold indices (see below)'
'- Sparse matrices are written dense.'
''
'## Indices are 0-based'
''
'MATLAB counts from 1. Every index of a unit, template, channel, community or seed written here'
'was converted to 0-based (value - 1) and is listed in the root attribute `index_valued`:'
''
'| function | index-valued data |'
'| --- | --- |'
'| CCG, computeConnectivityCCG | `inputs/spike_units` (unit of each spike); output `ccg` / `CCGs` axes 2 and 3 follow this unit order; `outputs/ExcitatoryConnection`, `outputs/InhibitoryConnection` ([source, target] unit pairs) |'
'| computeActivityFeatures, sttc (all-pairs cases) | `inputs/spike_units` |'
'| computeWaveformFeatures | `inputs/peak_channel` (channel of each unit''s template used for the waveform) |'
'| detectNetworkBursts | `inputs/spike_units`, `outputs/bursts/Superburst` (index into the superburst list, -1 = not part of a superburst), `outputs/info/Units` column `Unit` |'
'| propagateLabels | `inputs/seed_idx` (seed nodes) |'
'| communityEnrichmentTest | `inputs/M` (community of each unit), `outputs/T` column `Community`, `outputs/inh_comm`, `outputs/exc_comm` |'
'| zscoreByGroup | `outputs/fallback` (group indices) |'
''
'Not converted: labels without index meaning (`inputs/strata`, ComBat batch labels, `seed_y` classes'
'1 = excitatory / 2 = inhibitory), `Level` and `Scales` numbering in detectNetworkBursts tables,'
'and graph-theory outputs that count nodes (`N`, `K`, `Nk`, `n_iter`). Matrix rows / columns of'
'graphs are the node order; node k of a matrix is unit k - 1.'
''
'## Tables'
''
'MATLAB tables are written with `parquetwrite` to `<case>_<name>.parquet` (`<name>` = path below'
'`inputs/` or `outputs/`, `/` replaced by `_`) and listed in the root attribute `tables`. Index'
'columns are converted as above.'
''
'## Not exported'
''
'Stochastic functions are not exported: `community_louvain` (random node order), the surrogate based'
'significance of `computeConnectivitySTTC`, `randmio_*` rewiring and hence `computeGraphFeatures`'
'(Louvain modularity, small-worldness), and the bootstrap / permutation functions.'
    };
    fid = fopen(fullfile(outDir, 'README.md'), 'w');
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end

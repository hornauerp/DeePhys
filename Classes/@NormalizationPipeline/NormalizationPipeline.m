classdef NormalizationPipeline < handle
% NORMALIZATIONPIPELINE  Configurable, serializable normalization for ML pipelines.
%
%   Encapsulates multi-step normalization (per-group z-score, global z-score,
%   max-abs scaling). Fits on training data, applies the same transform to
%   test data — preventing data leakage.
%
%   USAGE:
%     np = NormalizationPipeline.groupThenGlobal();
%     [X_train, np] = np.fit_transform(X_train, group_labels_train);
%     X_test = np.transform(X_test, group_labels_test);
%
%   The FitParams property stores all fitted statistics, enabling
%   serialization and reapplication to new data (e.g., deployment).

    properties
        Steps    cell = {}   % Cell array of step structs: {struct('type','...','params',struct(...))}
        FitParams struct = struct()  % Stored after fit_transform: per-step means, stds, scale factors
    end

    methods

        function np = NormalizationPipeline(steps)
        % Constructor. steps is a cell array of step structs.
            arguments
                steps cell = {}
            end
            np.Steps = steps;
        end

        function [X_out, np] = fit_transform(np, X, group_labels, options)
        % FIT_TRANSFORM  Fit normalization on training data and return transformed result.
        %
        %   [X_out, np] = np.fit_transform(X_train, group_labels_train)
        %   [X_out, np] = np.fit_transform(X_train, group_labels_train, CovariateLabels=Y)
        %
        % INPUTS:
        %   X               - (N x F) numeric matrix
        %   group_labels    - (N x 1) group assignment vector (numeric or categorical).
        %                     Used by group_zscore and combat steps. Pass [] if neither
        %                     step is present.
        %   CovariateLabels - (N x 1) biological class labels of the TRAINING rows, used to
        %                     separate class composition from batch effects when ComBat
        %                     estimates them. Pass [] (default) to estimate without.
        %
        % OUTPUTS:
        %   X_out - normalized training matrix
        %   np    - pipeline with FitParams populated
            arguments
                np NormalizationPipeline
                X {isnumeric}
                group_labels = []
                options.CovariateLabels = []
            end

            np.FitParams = struct();
            X_out = X;

            for i = 1:length(np.Steps)
                step = np.Steps{i};
                switch step.type
                    case 'group_zscore'
                        [X_out, fit] = fitGroupZscore(X_out, group_labels);
                        np.FitParams(1).group_zscore = fit;

                    case 'global_zscore'
                        [X_out, gmu, gsd] = normalize(X_out);
                        X_out(:, gsd == 0) = 0;
                        np.FitParams(1).global_zscore.mean = gmu;
                        np.FitParams(1).global_zscore.sd = gsd;

                    case 'clip'
                        clip_range = step.params.range;
                        X_out = clipKeepNaN(X_out, clip_range);
                        np.FitParams(1).clip.range = clip_range;

                    case 'max_abs_scale'
                        sf = max(abs(X_out));
                        sf(sf == 0) = 1;
                        X_out = X_out ./ sf;
                        np.FitParams(1).max_abs_scale.scale_factor = sf;

                    case 'combat'
                        [X_out, fit] = fitComBat(X_out, group_labels, options.CovariateLabels);
                        np.FitParams(1).combat = fit;

                    otherwise
                        error('NormalizationPipeline:unknownStep', ...
                            'Unknown normalization step type: "%s".', step.type);
                end
            end
        end

        function X_out = transform(np, X, group_labels, options)
        % TRANSFORM  Apply previously fitted normalization to new data.
        %
        %   X_out = np.transform(X_test, group_labels_test)
        %
        % Never pass labels of the data being transformed: ComBat's correction is label-free
        % (labels only inform fitting), and CovariateLabels is ignored here. It is accepted
        % for backward compatibility only.
            arguments
                np NormalizationPipeline
                X {isnumeric}
                group_labels = []
                options.CovariateLabels = []
            end

            X_out = X;
            for i = 1:length(np.Steps)
                step = np.Steps{i};
                switch step.type
                    case 'group_zscore'
                        fit = np.FitParams.group_zscore;
                        X_out = applyGroupZscore(X_out, group_labels, fit);

                    case 'global_zscore'
                        gfit = np.FitParams.global_zscore;
                        X_out = normalize(X_out, 'center', gfit.mean, 'scale', gfit.sd);
                        X_out(:, gfit.sd == 0) = 0;

                    case 'clip'
                        clip_range = np.FitParams.clip.range;
                        X_out = clipKeepNaN(X_out, clip_range);

                    case 'max_abs_scale'
                        sf = np.FitParams.max_abs_scale.scale_factor;
                        X_out = X_out ./ sf;

                    case 'combat'
                        fit = np.FitParams.combat;
                        X_out = applyComBat(X_out, group_labels, options.CovariateLabels, fit);
                end
            end
        end

    end

    methods (Static)

        function np = groupThenGlobal()
        % GROUPTTHENGLOBAL  Per-group z-score → global z-score → clip → max-abs scale.
            np = NormalizationPipeline({
                struct('type', 'group_zscore',  'params', struct())
                struct('type', 'global_zscore', 'params', struct())
                struct('type', 'clip',          'params', struct('range', 5))
                struct('type', 'max_abs_scale', 'params', struct())
            });
        end

        function np = globalOnly()
        % GLOBALONLY  Global z-score → clip → max-abs scale (no per-group step).
            np = NormalizationPipeline({
                struct('type', 'global_zscore', 'params', struct())
                struct('type', 'clip',          'params', struct('range', 5))
                struct('type', 'max_abs_scale', 'params', struct())
            });
        end

        function np = combatThenGlobal()
        % COMBATTTHENGLOBAL  ComBat batch correction → global z-score → clip → max-abs scale.
        %
        %   Use when recordings come from multiple chips/batches with systematic
        %   technical offsets. Pass chip/batch assignments as group_labels and
        %   the biological variable (e.g., Mutation) as CovariateLabels to
        %   prevent biological signal from being removed.
            np = NormalizationPipeline({
                struct('type', 'combat',        'params', struct())
                struct('type', 'global_zscore', 'params', struct())
                struct('type', 'clip',          'params', struct('range', 5))
                struct('type', 'max_abs_scale', 'params', struct())
            });
        end

    end
end

% =========================================================================
% Local helper functions
% =========================================================================

function [X_out, fit] = fitGroupZscore(X, group_labels)
% Fit per-group z-score on training data. Returns transformed X and fit struct.
% Groups with < 2 samples cannot be z-scored; their means/sds are set to NaN
% so that transform() routes them to the pooled-statistics fallback (also applied here,
% so training and new rows of a singleton group are treated the same way).
% Labels are stored as strings: the statistics are keyed by group name, so new data from
% the same plate/chip maps back to them whatever type the caller used.
    group_labels = string(group_labels(:));
    groups = unique(group_labels);
    fit.groups = groups;
    fit.means = nan(length(groups), size(X, 2));
    fit.sds   = nan(length(groups), size(X, 2));
    X_out = X;
    singleton = false(size(X, 1), 1);
    for g = 1:length(groups)
        mask = (group_labels == groups(g));
        if sum(mask) < 2
            warning('NormalizationPipeline:singletonGroup', ...
                'Group has only %d sample(s) — group z-score skipped; pooled statistics used instead.', sum(mask));
            singleton = singleton | mask;
        else
            [X_out(mask,:), fit.means(g,:), fit.sds(g,:)] = normalize(X(mask,:));
            X_out(mask, fit.sds(g,:) == 0) = 0;     % constant in the group: no deviation (not NaN)
        end
    end
    valid_g = ~any(isnan(fit.means), 2);
    if any(singleton) && any(valid_g)
        fallback_sd = mean(fit.sds(valid_g, :), 1);
        fallback_sd(fallback_sd == 0) = 1;
        X_out(singleton, :) = normalize(X(singleton, :), 'center', mean(fit.means(valid_g, :), 1), 'scale', fallback_sd);
    end
end

function X_out = applyGroupZscore(X, group_labels, fit)
% Apply previously fitted per-group z-score to new data.
% Unseen groups and singleton groups (NaN means) fall back to pooled training statistics.
    group_labels = string(group_labels(:));
    fit.groups = string(fit.groups(:));
    X_out = X;
    known_mask = false(size(X, 1), 1);
    for g = 1:length(fit.groups)
        mask = (group_labels == fit.groups(g));
        if any(mask) && ~any(isnan(fit.means(g,:)))
            X_out(mask,:) = normalize(X(mask,:), 'center', fit.means(g,:), 'scale', fit.sds(g,:));
            X_out(mask, fit.sds(g,:) == 0) = 0;
            known_mask = known_mask | mask(:);
        end
        % Groups with NaN means (singletons during fit) remain in ~known_mask → pooled fallback
    end
    % Groups unseen during fit (e.g. the held-out recording in recording-grouped CV) are
    % z-scored with their own statistics. Per-group z-scoring is label-free, so this does not
    % leak; using training statistics instead would leave the group offset in the data.
    [new_groups, ~, gidx] = unique(group_labels(~known_mask));
    rows = find(~known_mask);
    for g = 1:numel(new_groups)
        r = rows(gidx == g);
        if numel(r) >= 2 && ~ismember(new_groups(g), fit.groups)
            [X_out(r, :), ~, sd_g] = normalize(X(r, :));
            X_out(r, sd_g == 0) = 0;
            known_mask(r) = true;
        end
    end
    % Singleton groups: use mean of valid fitted group statistics
    if any(~known_mask)
        warning('NormalizationPipeline:unseenGroup', ...
            '%d samples belong to unseen or singleton groups — using pooled training statistics.', sum(~known_mask));
        valid_g = ~any(isnan(fit.means), 2);
        if any(valid_g)
            fallback_mean = mean(fit.means(valid_g, :), 1);
            fallback_sd   = mean(fit.sds(valid_g, :), 1);
        else
            fallback_mean = zeros(1, size(X, 2));
            fallback_sd   = ones(1, size(X, 2));
        end
        fallback_sd(fallback_sd == 0) = 1;
        X_out(~known_mask,:) = normalize(X(~known_mask,:), 'center', fallback_mean, 'scale', fallback_sd);
    end
end

function [X_out, fit] = fitComBat(X, batch_labels, covariate_labels)
% Parametric empirical Bayes batch correction (Johnson, Li & Rabinovic, 2007).
%
% Model per feature g: x = alpha + C*beta + gamma_b + sqrt(delta_b) * sigma * e.
%   1. Fit batch and covariate effects jointly (OLS); sigma^2 is the within-batch
%      residual variance, alpha the size-weighted grand mean.
%   2. Standardise: s = (x - alpha - C*beta) / sigma.
%   3. Per batch, gamma_hat / delta_hat are the mean / variance of s. Their priors are
%      estimated across features within each batch (normal for gamma, inverse gamma
%      for delta) and the posteriors are solved iteratively.
%
% Covariate labels are used only to estimate the batch effects (so that a batch's class
% composition is not mistaken for a batch effect). The correction itself does not use
% labels — see applyComBat — so training and test rows are transformed identically and
% no class information enters the features.
%
% INPUTS:
%   X               - (N x F) data matrix
%   batch_labels    - (N x 1) batch assignments (numeric or categorical)
%   covariate_labels - (N x 1) biological variable to protect (may be [])

    [N, F] = size(X);
    batch_labels = string(batch_labels(:));   % numeric, categorical, cellstr or string
    batches = unique(batch_labels);
    n_batch = numel(batches);
    B = zeros(N, n_batch);
    for b = 1:n_batch
        B(:, b) = batch_labels == batches(b);
    end
    n_per_batch = sum(B, 1)';

    % Covariate dummies (no intercept: the batch columns span it)
    C = buildDesignMatrix(covariate_labels, N);
    C = C(:, 2:end);
    if ~isempty(C) && rank([B, C]) < n_batch + size(C, 2)
        warning('NormalizationPipeline:combatConfounded', ...
            ['ComBat covariate is confounded with batch (%d/%d independent columns); ' ...
             'it cannot be protected and is ignored.'], rank([B, C]), n_batch + size(C, 2));
        C = zeros(N, 0);
    end

    % Missing values stay missing: each feature is fitted on its observed rows, and the EB
    % priors (pooled across features) use only finite estimates, so one missing entry
    % cannot turn every feature NaN.
    D        = [B, C];
    obs      = ~isnan(X);
    beta_hat = zeros(size(D, 2), F);
    full_col = all(obs, 1);
    beta_hat(:, full_col) = D \ X(:, full_col);              % (n_batch + P) x F
    for f = find(~full_col)
        r = obs(:, f);
        if any(r)
            beta_hat(:, f) = lsqminnorm(D(r, :), X(r, f));   % batches without data get 0
        end
    end
    alpha    = (n_per_batch' / N) * beta_hat(1:n_batch, :);  % 1 x F grand mean
    beta_cov = beta_hat(n_batch+1:end, :);                   % P x F covariate effects
    var_pooled = mean((X - D * beta_hat).^2, 1, 'omitnan');  % within-batch variance
    var_pooled(var_pooled == 0 | ~isfinite(var_pooled)) = 1;
    sigma = sqrt(var_pooled);

    S = (X - alpha - C * beta_cov) ./ sigma;

    gamma_hat = zeros(n_batch, F);
    delta_hat = ones(n_batch, F);
    for b = 1:n_batch
        mask = B(:, b) > 0;
        gamma_hat(b, :) = mean(S(mask, :), 1, 'omitnan');
        delta_hat(b, :) = var(S(mask, :), 0, 1, 'omitnan');
    end
    gamma_hat(~isfinite(gamma_hat)) = 0;                     % no observed value in the batch
    delta_hat(delta_hat <= 0 | ~isfinite(delta_hat)) = 1;

    gamma_star = gamma_hat;
    delta_star = delta_hat;
    if n_batch < 2
        % Nothing to correct between batches; return the data unchanged.
        gamma_star(:) = 0;
        delta_star(:) = 1;
    else
        for b = 1:n_batch
            [gamma_star(b, :), delta_star(b, :)] = ebPosterior(S(B(:, b) > 0, :), ...
                gamma_hat(b, :), delta_hat(b, :));
        end
    end

    fit.batches    = batches;
    fit.alpha      = alpha;
    fit.sigma      = sigma;
    fit.gamma_star = gamma_star;
    fit.delta_star = delta_star;
    X_out = applyComBat(X, batch_labels, [], fit);
end

function X_out = applyComBat(X, batch_labels, ~, fit)
% Apply a fitted ComBat correction: x' = alpha + (x - alpha - sigma*gamma_b) / sqrt(delta_b).
% Identical for training and new data and label-free (the covariate argument is ignored;
% it is kept for call compatibility). Batches unseen during fit get their location and
% scale estimated from their own samples — transductive but label-free, which is the
% normal case for grouped cross-validation where the batch is a held-out recording.
    [N, F] = size(X);
    S = (X - fit.alpha) ./ fit.sigma;
    X_out = X;
    [new_batches, ~, idx] = unique(string(batch_labels(:)));
    for k = 1:numel(new_batches)
        mask = idx == k;
        b = find(string(fit.batches) == new_batches(k), 1);
        if ~isempty(b)
            g = fit.gamma_star(b, :);  d = fit.delta_star(b, :);
        elseif sum(mask) >= 2
            g = mean(S(mask, :), 1, 'omitnan');
            d = var(S(mask, :), 0, 1, 'omitnan');
            g(~isfinite(g)) = 0;
            d(d <= 0 | ~isfinite(d)) = 1;
            [g, d] = ebPosterior(S(mask, :), g, d);
        else
            warning('NormalizationPipeline:unseenBatch', ...
                'Batch "%s" has a single unseen sample — no batch correction applied.', string(new_batches(k)));
            continue
        end
        X_out(mask, :) = fit.alpha + (S(mask, :) - g) ./ sqrt(d) .* fit.sigma;
    end
    if F == 0, X_out = zeros(N, 0); end
end

function [g_star, d_star] = ebPosterior(S_b, g_hat, d_hat)
% Parametric EB posteriors for one batch (Johnson et al. 2007, it.sol): priors are
% estimated across features (gamma ~ N(g_bar, t2), delta ~ InvGamma(a, b)).
    n = sum(~isnan(S_b), 1);                 % observed samples per feature
    has = n >= 2;                            % features that inform the priors
    g_bar = mean(g_hat(has));
    t2    = var(g_hat(has));
    m     = mean(d_hat(has));
    s2    = var(d_hat(has));
    if sum(has) < 2 || t2 <= 0 || s2 <= 0 || size(S_b, 1) < 2
        g_star = g_hat;  d_star = d_hat;   % too little information to shrink
        return
    end
    a = (2 * s2 + m^2) / s2;
    b = (m * s2 + m^3) / s2;
    g_old = g_hat;  d_old = d_hat;
    for it = 1:100
        g_new = (n .* t2 .* g_hat + d_old * g_bar) ./ (n .* t2 + d_old);
        sum2  = sum((S_b - g_new).^2, 1, 'omitnan');
        d_new = (b + 0.5 * sum2) ./ (n / 2 + a - 1);
        change = max(max(abs(g_new - g_old) ./ max(abs(g_old), eps)), ...
                     max(abs(d_new - d_old) ./ max(abs(d_old), eps)));
        g_old = g_new;  d_old = d_new;
        if change < 1e-4, break; end
    end
    g_star = g_new;  d_star = max(d_new, 1e-6);
end

function Z = buildDesignMatrix(covariate_labels, N, cats)
% Build design matrix: intercept + one-hot coding of covariate_labels (drop first level).
% cats: cell of known categories (used at transform time to match training encoding).
    if isempty(covariate_labels)
        Z = ones(N, 1);
        return
    end
    if nargin < 3 || isempty(cats)
        cats = getCovariateCategories(covariate_labels);
    end
    n_cats = numel(cats);
    if n_cats <= 1
        Z = ones(N, 1);
        return
    end
    % Dummy coding: drop first category (reference level)
    Z = ones(N, n_cats);           % col 1 = intercept
    for k = 2:n_cats
        Z(:, k) = double(covariate_labels == cats{k});
    end
end

function cats = getCovariateCategories(covariate_labels)
% Return sorted unique categories for consistent dummy coding across fit/transform.
    if isempty(covariate_labels)
        cats = {};
        return
    end
    u = unique(covariate_labels);
    if iscell(u) || isstring(u)
        cats = cellstr(u);
    else
        cats = num2cell(u);
    end
end

function X = clipKeepNaN(X, r)
% Clip to [-r, r]. min/max would turn NaN (missing values, constant features) into +-r.
X(X > r)  = r;
X(X < -r) = -r;
end

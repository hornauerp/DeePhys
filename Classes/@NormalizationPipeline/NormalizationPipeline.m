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
        %   CovariateLabels - (N x 1) biological class labels to protect during ComBat
        %                     correction. Required when pipeline contains a 'combat' step.
        %                     Pass [] (default) to correct without protecting covariates.
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
                        np.FitParams(1).global_zscore.mean = gmu;
                        np.FitParams(1).global_zscore.sd = gsd;

                    case 'clip'
                        clip_range = step.params.range;
                        X_out = max(min(X_out, clip_range), -clip_range);
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
        %   X_out = np.transform(X_test, group_labels_test, CovariateLabels=Y_test)
        %
        % CovariateLabels is required when the pipeline contains a 'combat' step.
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

                    case 'clip'
                        clip_range = np.FitParams.clip.range;
                        X_out = max(min(X_out, clip_range), -clip_range);

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
% so that transform() routes them to the pooled-statistics fallback.
    groups = unique(group_labels);
    fit.groups = groups;
    fit.means = nan(length(groups), size(X, 2));
    fit.sds   = nan(length(groups), size(X, 2));
    X_out = X;
    for g = 1:length(groups)
        mask = (group_labels == groups(g));
        if sum(mask) < 2
            warning('NormalizationPipeline:singletonGroup', ...
                'Group has only %d sample(s) — group z-score skipped; pooled statistics used instead.', sum(mask));
        else
            [X_out(mask,:), fit.means(g,:), fit.sds(g,:)] = normalize(X(mask,:));
        end
    end
end

function X_out = applyGroupZscore(X, group_labels, fit)
% Apply previously fitted per-group z-score to new data.
% Unseen groups and singleton groups (NaN means) fall back to pooled training statistics.
    X_out = X;
    known_mask = false(size(X, 1), 1);
    for g = 1:length(fit.groups)
        mask = (group_labels == fit.groups(g));
        if any(mask) && ~any(isnan(fit.means(g,:)))
            X_out(mask,:) = normalize(X(mask,:), 'center', fit.means(g,:), 'scale', fit.sds(g,:));
            known_mask = known_mask | mask(:);
        end
        % Groups with NaN means (singletons during fit) remain in ~known_mask → pooled fallback
    end
    % Handle unseen groups and singleton groups: use mean of valid fitted group statistics
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
% Parametric empirical Bayes batch correction (Johnson et al., 2007).
%
% Estimates batch-specific additive (gamma) and multiplicative (delta) effects
% on standardized residuals after removing biological signal via OLS. Empirical
% Bayes shrinkage pools information across features to stabilize estimates for
% small batches.
%
% INPUTS:
%   X               - (N x F) data matrix
%   batch_labels    - (N x 1) batch assignments (numeric or categorical)
%   covariate_labels - (N x 1) biological variable to protect (may be [])

    [N, F] = size(X);
    batches = unique(batch_labels);
    n_batch = numel(batches);

    % Build design matrix for biological covariates (intercept + dummy coding)
    Z = buildDesignMatrix(covariate_labels, N);

    % OLS: estimate biological effects, then get residuals
    beta_hat = Z \ X;                     % (P x F)
    R = X - Z * beta_hat;                 % (N x F) batch-free residuals

    % Pooled per-feature variance across all samples
    var_pooled = var(R, 0, 1);            % (1 x F)
    var_pooled(var_pooled == 0 | isnan(var_pooled)) = 1;

    % Standardize residuals
    X_std = R ./ sqrt(var_pooled);        % (N x F)

    % Per-batch mean (gamma_hat) and variance (delta_hat) of standardized residuals
    gamma_hat = zeros(n_batch, F);
    delta_hat = ones(n_batch, F);
    n_per_batch = zeros(n_batch, 1);
    for b = 1:n_batch
        mask = (batch_labels == batches(b));
        n_per_batch(b) = sum(mask);
        if sum(mask) >= 1
            gamma_hat(b,:) = mean(X_std(mask,:), 1);
        end
        if sum(mask) >= 2
            delta_hat(b,:) = var(X_std(mask,:), 0, 1);
        end
    end
    delta_hat(delta_hat == 0 | isnan(delta_hat)) = 1;

    % Empirical Bayes hyperparameters (moment-matching across batches)
    gamma_bar = mean(gamma_hat, 1);       % grand mean per feature
    tau2_bar  = var(gamma_hat, 0, 1);    % between-batch variance per feature
    tau2_bar(tau2_bar == 0) = 1e-6;

    a_prior = (mean(delta_hat, 1).^2 + 2 .* var(delta_hat, 0, 1)) ./ ...
              max(var(delta_hat, 0, 1), 1e-10);   % shape
    b_prior = (mean(delta_hat, 1) .* var(delta_hat, 0, 1) + mean(delta_hat, 1).^3) ./ ...
              max(var(delta_hat, 0, 1), 1e-10);   % rate

    % Posterior estimates (shrinkage toward grand mean / prior)
    gamma_star = zeros(n_batch, F);
    delta_star = ones(n_batch, F);
    for b = 1:n_batch
        nb = n_per_batch(b);
        % Additive shrinkage: precision-weighted average of batch estimate and grand mean
        w = nb ./ (nb + var_pooled ./ max(tau2_bar, 1e-10));
        gamma_star(b,:) = w .* gamma_hat(b,:) + (1 - w) .* gamma_bar;

        % Multiplicative: inverse-gamma posterior mode
        delta_star(b,:) = (b_prior + nb/2 .* delta_hat(b,:)) ./ ...
                          (a_prior + nb/2 + 1);
    end
    delta_star = max(delta_star, 1e-6);   % keep positive

    % Apply correction to training data
    X_out = zeros(N, F);
    for b = 1:n_batch
        mask = (batch_labels == batches(b));
        if ~any(mask); continue; end
        X_corrected = (X_std(mask,:) - gamma_star(b,:)) ./ sqrt(delta_star(b,:));
        X_out(mask,:) = X_corrected .* sqrt(var_pooled) + Z(mask,:) * beta_hat;
    end

    % Store all parameters needed for applyComBat
    fit.batches            = batches;
    fit.beta_hat           = beta_hat;
    fit.var_pooled         = var_pooled;
    fit.gamma_star         = gamma_star;
    fit.delta_star         = delta_star;
    fit.covariate_cats     = getCovariateCategories(covariate_labels);
end

function X_out = applyComBat(X, batch_labels, covariate_labels, fit)
% Apply ComBat correction fitted on training data to new samples.
% Unseen batches receive no correction and trigger a warning.
    [N, ~] = size(X);

    % Reconstruct design matrix using training-time covariate categories
    Z = buildDesignMatrix(covariate_labels, N, fit.covariate_cats);

    % Recover biological signal estimate for new data
    X_hat = Z * fit.beta_hat;             % (N x F)
    R = X - X_hat;
    X_std = R ./ sqrt(fit.var_pooled);    % standardize using training variance

    X_out = X;
    seen = false(N, 1);
    for b = 1:numel(fit.batches)
        mask = (batch_labels == fit.batches(b));
        if ~any(mask); continue; end
        X_corrected = (X_std(mask,:) - fit.gamma_star(b,:)) ./ sqrt(fit.delta_star(b,:));
        X_out(mask,:) = X_corrected .* sqrt(fit.var_pooled) + X_hat(mask,:);
        seen = seen | mask(:);
    end

    unseen = ~seen;
    if any(unseen)
        warning('NormalizationPipeline:unseenBatch', ...
            '%d samples belong to batches unseen during fit — no batch correction applied.', sum(unseen));
    end
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

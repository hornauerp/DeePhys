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

        function [X_out, np] = fit_transform(np, X, group_labels)
        % FIT_TRANSFORM  Fit normalization on training data and return transformed result.
        %
        %   [X_out, np] = np.fit_transform(X_train, group_labels_train)
        %
        % INPUTS:
        %   X             - (N x F) numeric matrix
        %   group_labels  - (N x 1) group assignment vector (numeric or categorical).
        %                   Pass [] if no group normalization step is used.
        %
        % OUTPUTS:
        %   X_out - normalized training matrix
        %   np    - pipeline with FitParams populated
            arguments
                np NormalizationPipeline
                X {isnumeric}
                group_labels = []
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

                    otherwise
                        error('NormalizationPipeline:unknownStep', ...
                            'Unknown normalization step type: "%s".', step.type);
                end
            end
        end

        function X_out = transform(np, X, group_labels)
        % TRANSFORM  Apply previously fitted normalization to new data.
        %
        %   X_out = np.transform(X_test, group_labels_test)
            arguments
                np NormalizationPipeline
                X {isnumeric}
                group_labels = []
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

    end
end

% =========================================================================
% Local helper functions
% =========================================================================

function [X_out, fit] = fitGroupZscore(X, group_labels)
% Fit per-group z-score on training data. Returns transformed X and fit struct.
    groups = unique(group_labels);
    fit.groups = groups;
    fit.means = zeros(length(groups), size(X, 2));
    fit.sds   = zeros(length(groups), size(X, 2));
    X_out = X;
    for g = 1:length(groups)
        mask = (group_labels == groups(g));
        [X_out(mask,:), fit.means(g,:), fit.sds(g,:)] = normalize(X(mask,:));
    end
end

function X_out = applyGroupZscore(X, group_labels, fit)
% Apply previously fitted per-group z-score to new data.
% Unseen groups fall back to pooled training statistics.
    X_out = X;
    known_mask = false(size(X, 1), 1);
    for g = 1:length(fit.groups)
        mask = (group_labels == fit.groups(g));
        if any(mask)
            X_out(mask,:) = normalize(X(mask,:), 'center', fit.means(g,:), 'scale', fit.sds(g,:));
            known_mask = known_mask | mask(:);
        end
    end
    % Handle unseen groups: use mean of fitted group statistics as fallback
    if any(~known_mask)
        warning('NormalizationPipeline:unseenGroup', ...
            '%d samples belong to groups not seen during fit — using pooled training statistics.', sum(~known_mask));
        fallback_mean = mean(fit.means, 1);
        fallback_sd   = mean(fit.sds, 1);
        fallback_sd(fallback_sd == 0) = 1;
        X_out(~known_mask,:) = normalize(X(~known_mask,:), 'center', fallback_mean, 'scale', fallback_sd);
    end
end

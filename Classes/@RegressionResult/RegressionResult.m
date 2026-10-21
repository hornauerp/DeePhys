classdef RegressionResult < MLResult
% REGRESSIONRESULT  Result container for one cross-validation fold of regression.
%
%   Stores regression model, predictions, feature importance, and all
%   parameters used — decoupled from RecordingGroup state.
%
%   Inherits shared properties (Mdl, Y_pred, Y_test, objects, train_acc,
%   predImp, Features, Parameters) and summarizeImportance from MLResult.

    properties
        mse_train   % Training mean squared error (resubstitution loss)
    end

    methods

        function r = RegressionResult(varargin)
        % Constructor. Accepts name-value pairs or no arguments.
            r = r@MLResult(varargin{:});   % parent handles shared properties
            if nargin > 0
                for i = 1:2:length(varargin)
                    if strcmp(varargin{i}, 'mse_train')
                        r.mse_train = varargin{i+1};
                    end
                end
            end
        end

        function metrics = computeMetrics(r)
        % COMPUTEMETRICS  Compute MSE, MAE, R², correlation from predictions.
            Y_pred    = r.Y_pred(:);
            Y_test    = r.Y_test(:);
            residuals = Y_pred - Y_test;

            metrics.MSE  = mean(residuals.^2);
            metrics.RMSE = sqrt(metrics.MSE);
            metrics.MAE  = mean(abs(residuals));

            ss_tot = sum((Y_test - mean(Y_test)).^2);
            if ss_tot == 0
                metrics.R2 = NaN;
            else
                metrics.R2 = 1 - sum(residuals.^2) / ss_tot;
            end
            metrics.Correlation = corr(Y_pred, Y_test);
            metrics.N           = length(Y_test);
            metrics.mse_train   = r.mse_train;
        end

    end

    methods (Static)

        function T = results2table(result_array)
        % RESULTS2TABLE  Convert a RegressionResult array to a flat table.
        %
        % Returns a table with one row per test object across all folds:
        %   Fold, Y_pred, Y_test, Residual
            rows = {};
            for k = 1:length(result_array)
                r        = result_array(k);
                n        = length(r.Y_test);
                fold_col = repmat(k, n, 1);
                residual = r.Y_pred(:) - r.Y_test(:);

                rows{k} = table(fold_col, r.Y_pred(:), r.Y_test(:), residual, ...
                    'VariableNames', {'Fold', 'Y_pred', 'Y_test', 'Residual'}); %#ok<AGROW>
            end
            T = vertcat(rows{:});
        end

        function summary = summarizeFolds(result_array)
        % SUMMARIZEFOLDS  Compute per-fold and aggregate regression metrics.
        %
        % Returns struct with: per_fold (table), mean_MSE, mean_R2, mean_Correlation, and
        % R2_pooled: 1 - SS_res / SS_tot over the pooled out-of-fold predictions. The fold R²
        % compares each fold with its own mean; with grouped folds a fold can hold a single
        % target value (R² undefined) and few rows (unstable), so R2_pooled is the summary
        % to report and to test (Regressor.permutationTest). mean_R2 skips undefined folds.
            n_folds = length(result_array);
            mses = zeros(n_folds, 1);
            r2s  = zeros(n_folds, 1);
            cors = zeros(n_folds, 1);
            for k = 1:n_folds
                m      = result_array(k).computeMetrics();
                mses(k) = m.MSE;
                r2s(k)  = m.R2;
                cors(k) = m.Correlation;
            end
            summary.per_fold         = table((1:n_folds)', mses, r2s, cors, ...
                'VariableNames', {'Fold', 'MSE', 'R2', 'Correlation'});
            summary.mean_MSE         = mean(mses);
            summary.mean_R2          = mean(r2s, 'omitnan');
            summary.mean_Correlation = mean(cors, 'omitnan');
            y_true = cell2mat(arrayfun(@(r) r.Y_test(:), result_array(:), 'UniformOutput', false));
            y_pred = cell2mat(arrayfun(@(r) r.Y_pred(:), result_array(:), 'UniformOutput', false));
            ss_tot = sum((y_true - mean(y_true)).^2);
            if isempty(y_true) || ss_tot == 0
                summary.R2_pooled = NaN;
            else
                summary.R2_pooled = 1 - sum((y_pred - y_true).^2) / ss_tot;
            end
        end

    end
end

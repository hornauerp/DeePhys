classdef RegressionResult
% REGRESSIONRESULT  Standalone container for regression outputs.
%
%   Stores regression model, predictions, feature importance, and all
%   parameters used — decoupled from RecordingGroup state.

    properties
        Mdl             % Trained regression model
        Y_pred          % Predicted values (test set)
        Y_test          % True values (test set)
        objects         % Test objects (Unit/MEArecording/Culture array)
        mse_train       % Training mean squared error
        train_acc       % Training accuracy (R² or correlation)
        predImp         % Feature importance
        Features        % Feature names used (string array)
        Parameters      % Struct of all parameters: algorithm, normalization, CV, etc.
    end

    methods
        function r = RegressionResult(varargin)
            if nargin > 0
                for i = 1:2:length(varargin)
                    r.(varargin{i}) = varargin{i+1};
                end
            end
        end

        function metrics = computeMetrics(r)
        % COMPUTEMETRICS  Compute MSE, MAE, R², correlation from predictions.
            Y_pred = r.Y_pred(:);
            Y_test = r.Y_test(:);
            residuals = Y_pred - Y_test;
            metrics.MSE = mean(residuals.^2);
            metrics.RMSE = sqrt(metrics.MSE);
            metrics.MAE = mean(abs(residuals));
            ss_tot = sum((Y_test - mean(Y_test)).^2);
            if ss_tot == 0
                metrics.R2 = NaN;   % undefined when all Y_test values are identical
            else
                metrics.R2 = 1 - sum(residuals.^2) / ss_tot;
            end
            metrics.Correlation = corr(Y_pred, Y_test);
            metrics.N = length(Y_test);
            metrics.mse_train = r.mse_train;
        end
    end

    methods (Static)

        function T = results2table(result_array)
        % RESULTS2TABLE  Convert an array of RegressionResults to a flat table.
        %
        %   T = RegressionResult.results2table(result_array)
        %
        % Returns a table with one row per test object across all folds.
            rows = {};
            for k = 1:length(result_array)
                r = result_array(k);
                n = length(r.Y_test);
                fold_col = repmat(k, n, 1);
                residual = r.Y_pred(:) - r.Y_test(:);

                fold_table = table(fold_col, r.Y_pred(:), r.Y_test(:), residual, ...
                    'VariableNames', {'Fold', 'Y_pred', 'Y_test', 'Residual'});
                rows{k} = fold_table;
            end
            T = vertcat(rows{:});
        end

        function summary = summarizeFolds(result_array)
        % SUMMARIZEFOLDS  Compute per-fold and aggregate regression metrics.
            n_folds = length(result_array);
            mses = zeros(n_folds, 1);
            r2s  = zeros(n_folds, 1);
            cors = zeros(n_folds, 1);
            for k = 1:n_folds
                m = result_array(k).computeMetrics();
                mses(k) = m.MSE;
                r2s(k)  = m.R2;
                cors(k) = m.Correlation;
            end
            summary.per_fold = table((1:n_folds)', mses, r2s, cors, ...
                'VariableNames', {'Fold', 'MSE', 'R2', 'Correlation'});
            summary.mean_MSE = mean(mses);
            summary.mean_R2 = mean(r2s);
            summary.mean_Correlation = mean(cors);
        end

    end
end

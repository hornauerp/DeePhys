classdef ClassificationResult < MLResult
% CLASSIFICATIONRESULT  Result container for one cross-validation fold of classification.
%
%   Stores classifier model, predictions, scores, feature importance, and
%   all parameters used — decoupled from RecordingGroup state.
%
%   Inherits shared properties (Mdl, Y_pred, Y_test, objects, train_acc,
%   predImp, Features, Parameters) and summarizeImportance from MLResult.

    properties
        scores      % Prediction scores / posterior probabilities
        GroupLabels % Class label names used in this fold
    end

    methods

        function r = ClassificationResult(varargin)
        % Constructor. Accepts name-value pairs or no arguments.
            r = r@MLResult(varargin{:});   % parent handles shared properties
            if nargin > 0
                child_props = {'scores', 'GroupLabels'};
                for i = 1:2:length(varargin)
                    if ismember(varargin{i}, child_props)
                        r.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end

        function metrics = computeMetrics(r)
        % COMPUTEMETRICS  Accuracy and macro-averaged precision, recall and F1 (over the
        %   classes present in Y_test or Y_pred). Always macro: the F1 of one arbitrary
        %   "positive" class for two-class folds made mean_F1 average different metrics when
        %   grouped folds of a multiclass problem hold two classes only.
            Y_pred = r.Y_pred(:);
            Y_test = r.Y_test(:);
            metrics.Accuracy  = mean(Y_pred == Y_test);
            metrics.N         = length(Y_test);
            metrics.train_acc = r.train_acc;
            [metrics.Precision, metrics.Recall, metrics.F1_score] = ...
                ClassificationResult.macroScores(Y_pred, Y_test);
            metrics.ConfusionMatrix = confusionmat(Y_test, Y_pred);
        end

    end

    methods (Static)

        function [P, R, F1] = macroScores(Y_pred, Y_test)
        % MACROSCORES  Macro-averaged precision, recall and F1 over the classes in Y_test ∪ Y_pred.
            classes = unique([Y_pred(:); Y_test(:)]);
            n = numel(classes);
            Pc = zeros(n, 1); Rc = Pc; Fc = Pc;
            for c = 1:n
                tp = sum(Y_pred == classes(c) & Y_test == classes(c));
                fp = sum(Y_pred == classes(c) & Y_test ~= classes(c));
                fn = sum(Y_pred ~= classes(c) & Y_test == classes(c));
                Pc(c) = tp / max(tp + fp, 1);
                Rc(c) = tp / max(tp + fn, 1);
                Fc(c) = 2 * Pc(c) * Rc(c) / max(Pc(c) + Rc(c), eps);
            end
            P = mean(Pc); R = mean(Rc); F1 = mean(Fc);
        end

        function T = results2table(result_array)
        % RESULTS2TABLE  Convert a ClassificationResult array to a flat table.
        %
        % Returns a table with one row per test object across all folds:
        %   Fold, Y_pred, Y_test, Correct, Score_1, Score_2, ...
            rows = {};
            for k = 1:length(result_array)
                r = result_array(k);
                n = length(r.Y_test);
                fold_col = repmat(k, n, 1);
                correct  = r.Y_pred(:) == r.Y_test(:);

                fold_table = table(fold_col, r.Y_pred(:), r.Y_test(:), correct, ...
                    'VariableNames', {'Fold', 'Y_pred', 'Y_test', 'Correct'});

                if ~isempty(r.scores)
                    score_names = "Score_" + (1:size(r.scores, 2));
                    score_table = array2table(r.scores, 'VariableNames', score_names);
                    fold_table  = [fold_table, score_table]; %#ok<AGROW>
                end
                rows{k} = fold_table; %#ok<AGROW>
            end
            T = vertcat(rows{:});
        end

        function summary = summarizeFolds(result_array)
        % SUMMARIZEFOLDS  Compute per-fold and aggregate classification metrics.
        %
        % Returns struct with: per_fold (table), mean_accuracy, std_accuracy, mean_F1, and
        % balanced_accuracy (mean per-class recall over the pooled out-of-fold predictions,
        % robust to class imbalance and to unequal grouped-CV fold sizes) and macro_F1 (macro-F1
        % over the same pooled predictions).
            n_folds = length(result_array);
            accs = zeros(n_folds, 1);
            f1s  = zeros(n_folds, 1);
            for k = 1:n_folds
                m      = result_array(k).computeMetrics();
                accs(k) = m.Accuracy;
                f1s(k)  = m.F1_score;
            end
            summary.per_fold      = table((1:n_folds)', accs, f1s, ...
                'VariableNames', {'Fold', 'Accuracy', 'F1_score'});
            summary.mean_accuracy = mean(accs);
            summary.std_accuracy  = std(accs);
            summary.mean_F1       = mean(f1s);
            summary.std_F1        = std(f1s);
            y_true = cell2mat(arrayfun(@(r) r.Y_test(:), result_array(:), 'UniformOutput', false));
            y_pred = cell2mat(arrayfun(@(r) r.Y_pred(:), result_array(:), 'UniformOutput', false));
            classes = unique(y_true);
            recall = arrayfun(@(c) mean(y_pred(y_true == c) == c), classes);
            summary.balanced_accuracy = mean(recall);
            [~, ~, summary.macro_F1] = ClassificationResult.macroScores(y_pred, y_true);   % pooled out-of-fold
        end

    end
end

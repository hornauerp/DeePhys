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
        % COMPUTEMETRICS  Compute accuracy, F1, precision, recall from predictions.
            Y_pred = r.Y_pred(:);
            Y_test = r.Y_test(:);
            classes   = unique([Y_pred; Y_test]);
            n_classes = length(classes);

            metrics.Accuracy  = mean(Y_pred == Y_test);
            metrics.N         = length(Y_test);
            metrics.train_acc = r.train_acc;

            if n_classes == 2
                pos = classes(1);
                tp  = sum(Y_pred == pos & Y_test == pos);
                fp  = sum(Y_pred == pos & Y_test ~= pos);
                fn  = sum(Y_pred ~= pos & Y_test == pos);
                precision = tp / max(tp + fp, 1);
                recall    = tp / max(tp + fn, 1);
                metrics.Precision = precision;
                metrics.Recall    = recall;
                metrics.F1_score  = 2 * precision * recall / max(precision + recall, eps);
            else
                prec_sum = 0; rec_sum = 0; f1_sum = 0;
                for c = 1:n_classes
                    tp  = sum(Y_pred == classes(c) & Y_test == classes(c));
                    fp  = sum(Y_pred == classes(c) & Y_test ~= classes(c));
                    fn  = sum(Y_pred ~= classes(c) & Y_test == classes(c));
                    P_c = tp / max(tp + fp, 1);
                    R_c = tp / max(tp + fn, 1);
                    prec_sum = prec_sum + P_c;
                    rec_sum  = rec_sum  + R_c;
                    f1_sum   = f1_sum   + 2 * P_c * R_c / max(P_c + R_c, eps);
                end
                metrics.Precision = prec_sum / n_classes;
                metrics.Recall    = rec_sum  / n_classes;
                metrics.F1_score  = f1_sum   / n_classes;
            end
            metrics.ConfusionMatrix = confusionmat(Y_test, Y_pred);
        end

    end

    methods (Static)

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
        % Returns struct with: per_fold (table), mean_accuracy, std_accuracy, mean_F1.
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
        end

    end
end

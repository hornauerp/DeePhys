classdef ClassificationResult
% CLASSIFICATIONRESULT  Standalone container for classification outputs.
%
%   Stores classifier model, predictions, scores, feature importance, and
%   all parameters used — decoupled from RecordingGroup state.
%
%   Field names match the existing struct-based result format for backward
%   compatibility with scripts that access result.Y_pred, result.scores, etc.

    properties
        Mdl             % Trained classifier model (compact or full)
        Y_pred          % Predicted labels (test set)
        Y_test          % True labels (test set)
        scores          % Prediction scores / posterior probabilities
        objects         % Test objects (Unit/MEArecording/Culture array)
        train_acc       % Training accuracy
        predImp         % Feature importance (predictor importance table)
        GroupLabels     % Group/class labels used for classification
        Features        % Feature names used (string array)
        Parameters      % Struct of all parameters: algorithm, normalization, CV, etc.
    end

    methods
        function r = ClassificationResult(varargin)
        % Constructor. Accepts name-value pairs or no arguments.
            if nargin > 0
                for i = 1:2:length(varargin)
                    r.(varargin{i}) = varargin{i+1};
                end
            end
        end

        function metrics = computeMetrics(r)
        % COMPUTEMETRICS  Compute accuracy, F1, precision, recall from predictions.
            Y_pred = r.Y_pred(:);
            Y_test = r.Y_test(:);
            classes = unique([Y_pred; Y_test]);
            n_classes = length(classes);

            metrics.Accuracy = mean(Y_pred == Y_test);
            metrics.N = length(Y_test);
            metrics.train_acc = r.train_acc;

            if n_classes == 2
                % Binary: compute for the positive class (first label)
                pos = classes(1);
                tp = sum(Y_pred == pos & Y_test == pos);
                fp = sum(Y_pred == pos & Y_test ~= pos);
                fn = sum(Y_pred ~= pos & Y_test == pos);
                precision = tp / max(tp + fp, 1);
                recall    = tp / max(tp + fn, 1);
                metrics.Precision = precision;
                metrics.Recall = recall;
                metrics.F1_score = 2 * precision * recall / max(precision + recall, eps);
            else
                % Multiclass: macro-averaged
                prec_sum = 0; rec_sum = 0;
                for c = 1:n_classes
                    tp = sum(Y_pred == classes(c) & Y_test == classes(c));
                    fp = sum(Y_pred == classes(c) & Y_test ~= classes(c));
                    fn = sum(Y_pred ~= classes(c) & Y_test == classes(c));
                    prec_sum = prec_sum + tp / max(tp + fp, 1);
                    rec_sum  = rec_sum  + tp / max(tp + fn, 1);
                end
                metrics.Precision = prec_sum / n_classes;
                metrics.Recall = rec_sum / n_classes;
                metrics.F1_score = 2 * metrics.Precision * metrics.Recall / max(metrics.Precision + metrics.Recall, eps);
            end
            metrics.ConfusionMatrix = confusionmat(Y_test, Y_pred);
        end
    end

    methods (Static)

        function T = results2table(result_array)
        % RESULTS2TABLE  Convert an array of ClassificationResults to a flat table.
        %
        %   T = ClassificationResult.results2table(result_array)
        %
        % Returns a table with one row per test object across all folds:
        %   Fold, Y_pred, Y_test, Correct, Score_1, Score_2, ...
        %
        % Also appends fold-level metrics (Accuracy, F1) and parameters.
            rows = {};
            for k = 1:length(result_array)
                r = result_array(k);
                n = length(r.Y_test);
                fold_col = repmat(k, n, 1);
                correct = r.Y_pred(:) == r.Y_test(:);

                fold_table = table(fold_col, r.Y_pred(:), r.Y_test(:), correct, ...
                    'VariableNames', {'Fold', 'Y_pred', 'Y_test', 'Correct'});

                % Add scores if available
                if ~isempty(r.scores)
                    score_names = "Score_" + (1:size(r.scores, 2));
                    score_table = array2table(r.scores, 'VariableNames', score_names);
                    fold_table = [fold_table, score_table]; %#ok<AGROW>
                end

                rows{k} = fold_table;
            end
            T = vertcat(rows{:});
        end

        function summary = summarizeFolds(result_array)
        % SUMMARIZEFOLDS  Compute per-fold and aggregate metrics.
        %
        %   summary = ClassificationResult.summarizeFolds(result_array)
        %
        % Returns struct with: per_fold (table), mean_accuracy, mean_F1, std_accuracy.
            n_folds = length(result_array);
            accs = zeros(n_folds, 1);
            f1s  = zeros(n_folds, 1);
            for k = 1:n_folds
                m = result_array(k).computeMetrics();
                accs(k) = m.Accuracy;
                f1s(k)  = m.F1_score;
            end
            summary.per_fold = table((1:n_folds)', accs, f1s, ...
                'VariableNames', {'Fold', 'Accuracy', 'F1_score'});
            summary.mean_accuracy = mean(accs);
            summary.std_accuracy  = std(accs);
            summary.mean_F1       = mean(f1s);
            summary.std_F1        = std(f1s);
        end

    end
end

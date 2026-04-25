classdef MLResult
% MLRESULT  Shared base class for ML pipeline result containers.
%
% Inherited by ClassificationResult and RegressionResult.
% Holds the 8 properties common to both result types and the
% cross-class summarizeImportance utility.
%
% Do not instantiate directly — use ClassificationResult or RegressionResult.

    properties
        Mdl         % Trained model (compact or full)
        Y_pred      % Predicted values / labels (test set)
        Y_test      % True values / labels (test set)
        objects     % Test objects (Unit/MEArecording/Culture array)
        train_acc   % Training accuracy (classification) or R² (regression)
        predImp     % OOB predictor importance (RF only; [] otherwise)
        Features    % Feature names used (string array)
        Parameters  % Struct: algorithm, normalization_var, K_fold, N_hyper, seed
    end

    methods

        function r = MLResult(varargin)
        % Constructor — accepts name-value pairs for shared properties.
        % Called via r@MLResult(varargin{:}) from subclass constructors.
            if nargin > 0
                shared = {'Mdl','Y_pred','Y_test','objects','train_acc', ...
                          'predImp','Features','Parameters'};
                for i = 1:2:length(varargin)
                    if ismember(varargin{i}, shared)
                        r.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end

    end

    methods (Static)

        function T = summarizeImportance(result_array)
        % SUMMARIZEIMPORTANCE  Aggregate OOB predictor importance across CV folds.
        %
        % Works for both ClassificationResult and RegressionResult arrays.
        % Returns a table sorted by mean importance (descending):
        %   Feature | MeanImportance | StdImportance | MinImportance | MaxImportance
        %
        % Returns empty table when no fold has importance data (e.g., non-RF algorithm).
        %
        % EXAMPLE:
        %   result = exp.classify('Unit', 'Mutation', opts);
        %   T = ClassificationResult.summarizeImportance(result);
        %   disp(T(1:10, :));   % top 10 features
            valid_folds = find(~cellfun(@isempty, {result_array.predImp}));
            if isempty(valid_folds)
                T = table();
                return
            end

            all_names = string.empty;
            for k = valid_folds
                all_names = union(all_names, string(result_array(k).Features), 'stable');
            end
            n_feat  = numel(all_names);
            imp_mat = NaN(numel(valid_folds), n_feat);

            for fi = 1:numel(valid_folds)
                k     = valid_folds(fi);
                names = string(result_array(k).Features);
                imp   = result_array(k).predImp(:)';
                [~, loc] = ismember(names, all_names);
                imp_mat(fi, loc) = imp;
            end

            mean_imp = mean(imp_mat, 1, 'omitnan');
            std_imp  = std(imp_mat,  0, 1, 'omitnan');
            min_imp  = min(imp_mat,  [], 1);
            max_imp  = max(imp_mat,  [], 1);

            T = table(all_names(:), mean_imp(:), std_imp(:), min_imp(:), max_imp(:), ...
                'VariableNames', {'Feature','MeanImportance','StdImportance', ...
                                  'MinImportance','MaxImportance'});
            T = sortrows(T, 'MeanImportance', 'descend');
        end

    end
end

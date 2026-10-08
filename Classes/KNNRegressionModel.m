classdef KNNRegressionModel
% KNNREGRESSIONMODEL  Distance-weighted k-nearest-neighbour regression.
%
% MATLAB has fitcknn (classification) but no k-NN regression; MLPipeline.createRegressor
% returns this model for Algorithm = "knn". It supports the calls Regressor makes:
% predict(mdl, X) and resubLoss(mdl).
%
%   mdl  = KNNRegressionModel(X_train, Y_train, K);
%   yhat = predict(mdl, X_test);
%   K    = KNNRegressionModel.selectK(X_train, Y_train)   % leave-one-out choice of K

    properties
        X (:,:) double     % training features (already normalised upstream)
        Y (:,1) double     % training targets
        K (1,1) double = 5
    end

    methods
        function m = KNNRegressionModel(X, Y, K)
            if nargin < 3, K = 5; end
            m.X = X; m.Y = Y(:); m.K = min(K, size(X, 1));
        end

        function yhat = predict(m, Xq)
        % Inverse-distance weighted mean of the K nearest training targets.
            [idx, d] = knnsearch(m.X, Xq, 'K', m.K);
            w = 1 ./ max(d, 1e-12);
            yhat = sum(w .* reshape(m.Y(idx), size(idx)), 2) ./ sum(w, 2);
        end

        function L = resubLoss(m)
        % Mean squared error on the training data (optimistic, as for other models).
            L = mean((predict(m, m.X) - m.Y) .^ 2);
        end
    end

    methods (Static)
        function K = selectK(X, Y, candidates)
        % Leave-one-out mean squared error over candidate K values.
            if nargin < 3, candidates = [3 5 7 10 15 20 30]; end
            Y = Y(:);
            candidates = candidates(candidates < size(X, 1));
            if isempty(candidates), K = 1; return; end
            [idx, d] = knnsearch(X, X, 'K', max(candidates) + 1);
            idx = idx(:, 2:end); d = d(:, 2:end);            % drop the point itself
            w = 1 ./ max(d, 1e-12);
            Yn = reshape(Y(idx), size(idx));
            mse = arrayfun(@(k) mean((sum(w(:, 1:k) .* Yn(:, 1:k), 2) ./ sum(w(:, 1:k), 2) - Y) .^ 2), ...
                candidates);
            [~, best] = min(mse);
            K = candidates(best);
        end
    end
end

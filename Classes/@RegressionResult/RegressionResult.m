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
    end
end

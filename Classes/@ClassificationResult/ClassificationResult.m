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
    end
end

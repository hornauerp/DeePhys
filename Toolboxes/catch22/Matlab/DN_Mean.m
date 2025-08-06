function out = DN_Mean(y)

    % no combination of single functions
    coder.inline('never');
    
    out = mean(y);
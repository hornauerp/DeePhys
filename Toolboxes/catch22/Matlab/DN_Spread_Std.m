function out = DN_Spread_Std(y)

    % no combination of single functions
    coder.inline('never');
    
    out = std(y);
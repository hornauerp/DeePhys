function dt = parseDateRobust(date_val)
% PARSEDATEROBUST  Parse a date string in yyMMdd or yyyyMMdd format.
%
%   dt = parseDateRobust('221207')     % → 2022-12-07
%   dt = parseDateRobust('20221207')   % → 2022-12-07
%   dt = parseDateRobust(221207)       % numeric input converted to string
%
%   Returns NaT with a warning if the format cannot be determined.

    date_str = strtrim(string(date_val));

    if strlength(date_str) == 8
        dt = datetime(date_str, "InputFormat", "yyyyMMdd");
    elseif strlength(date_str) == 6
        dt = datetime(date_str, "InputFormat", "yyMMdd");
    else
        warning('DeePhys:UnknownDateFormat', ...
            'Cannot parse date "%s" — expected 6 (yyMMdd) or 8 (yyyyMMdd) digits.', date_str);
        dt = NaT;
    end
end

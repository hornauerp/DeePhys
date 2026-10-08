function dt = parseDateRobust(date_val)
% PARSEDATEROBUST  Parse a date given as yyMMdd, yyyyMMdd, yyyy-MM-dd or a datetime.
%
%   dt = parseDateRobust('221207')       % → 2022-12-07
%   dt = parseDateRobust('20221207')     % → 2022-12-07
%   dt = parseDateRobust('2022-12-07')   % → 2022-12-07 ('/' or '.' separators too)
%   dt = parseDateRobust(221207)         % numeric input converted to string
%
%   Returns NaT with a warning if the format cannot be determined.

    if isdatetime(date_val)
        dt = date_val;
        return
    end
    date_str = strtrim(string(date_val));
    dt = NaT;

    iso = regexp(date_str, '^(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})$', 'tokens', 'once');
    if ~isempty(iso)
        ymd = double(string(iso));
    elseif ~isempty(regexp(date_str, '^\d{8}$', 'once'))
        ymd = double([extractBetween(date_str, 1, 4), extractBetween(date_str, 5, 6), ...
                      extractBetween(date_str, 7, 8)]);
    elseif ~isempty(regexp(date_str, '^\d{6}$', 'once'))
        ymd = double([extractBetween(date_str, 1, 2), extractBetween(date_str, 3, 4), ...
                      extractBetween(date_str, 5, 6)]);
        % two-digit year: MATLAB's default pivot (the 100 years from now - 50)
        pivot = year(datetime('now')) - 50;
        ymd(1) = ymd(1) + 100 * floor((pivot - ymd(1)) / 100);
        if ymd(1) < pivot, ymd(1) = ymd(1) + 100; end
    else
        ymd = [];
    end

    if ~isempty(ymd) && ymd(2) >= 1 && ymd(2) <= 12 && ymd(3) >= 1 && ...
            ymd(3) <= eomday(ymd(1), ymd(2))
        dt = datetime(ymd(1), ymd(2), ymd(3));
    else
        warning('DeePhys:UnknownDateFormat', ...
            'Cannot parse date "%s": expected yyMMdd, yyyyMMdd or yyyy-MM-dd.', date_str);
    end
end

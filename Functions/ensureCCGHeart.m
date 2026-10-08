function ensureCCGHeart()
% ENSURECCGHEART  Compile the CCGHeart MEX used by CCG if this platform has no binary.
%
% The repository ships CCGHeart.c and a Windows binary (.mexw64); on other platforms
% the MEX is built on first use, which needs a C compiler ("mex -setup C").
    persistent ready
    if ~isempty(ready)
        return
    end
    func_dir = fileparts(mfilename('fullpath'));
    % Check the binary itself: which('CCGHeart') would also find the .c source.
    if ~isfile(fullfile(func_dir, ['CCGHeart.' mexext()]))
        prev_dir = cd(func_dir);
        restore  = onCleanup(@() cd(prev_dir));
        try
            mex('CCGHeart.c');
        catch ME
            error('DeePhys:CCGHeartMissing', ['CCGHeart.%s is missing and could not be ' ...
                'compiled (%s). Configure a C compiler with "mex -setup C" and retry.'], ...
                mexext(), ME.message);
        end
        clear restore
    end
    ready = true;
end

classdef CacheManager < handle
% CACHEMANAGER  Session-scoped disk cache for expensive computations.
%
% Singleton handle class. Manages a temp directory that persists for the
% lifetime of the MATLAB session and is automatically cleaned up when
% MATLAB exits (via the handle class destructor).
%
% Primarily used by CellTypeClassifier.getOrExtract to cache waveform/ACG
% extraction results on disk, so recreating a CellTypeClassifier with the
% same data and parameters avoids re-extraction.
%
% The directory is unique per MATLAB process (via PID), so concurrent MATLAB
% sessions do not interfere.
%
% USAGE:
%   cm   = CacheManager.instance();
%   path = cm.get('some_hash_key');   % returns full path to cache .mat file
%   if exist(path, 'file')
%       s = load(path, 'wf', 'acg', 'sr');
%   else
%       % compute...
%       save(path, 'wf', 'acg', 'sr', '-v7.3');
%   end
%
%   CacheManager.instance().CacheDir   % inspect cache directory
%   delete(CacheManager.instance())    % manually clean up (rarely needed)

    properties (SetAccess = private)
        CacheDir    string  % Absolute path to session temp directory
    end

    methods (Access = private)

        function cm = CacheManager()
            pid = feature('getpid');
            cm.CacheDir = fullfile(tempdir, sprintf('DeePhys_cache_%d', pid));
            if ~exist(cm.CacheDir, 'dir')
                mkdir(cm.CacheDir);
            end
        end

    end

    methods

        function path = get(cm, key)
        % GET  Return full path for a cache entry keyed by the given hash string.
            path = fullfile(cm.CacheDir, key + ".mat");
        end

        function delete(cm)
        % DELETE  Destructor — removes the cache directory.
        % Called automatically when MATLAB exits or when explicitly invoked.
            if exist(cm.CacheDir, 'dir')
                try
                    rmdir(cm.CacheDir, 's');
                catch
                end
            end
        end

    end

    methods (Static)

        function cm = instance()
        % INSTANCE  Return the singleton CacheManager for this MATLAB session.
            persistent shared
            if isempty(shared) || ~isvalid(shared)
                shared = CacheManager();
            end
            cm = shared;
        end

    end
end

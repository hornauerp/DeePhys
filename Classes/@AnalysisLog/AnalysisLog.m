classdef AnalysisLog < handle
% ANALYSISLOG  Lightweight provenance tracking for analysis sessions.
%
%   Records what analyses were run, with what parameters, and a summary of
%   results. Entries are appended in chronological order and can be exported
%   to a table or saved to a JSON file.
%
%   USAGE:
%     log = AnalysisLog();
%     log.add('classifyByFeatureGroups', struct('alg','rf','K_fold',5), 'Accuracy: 0.85');
%     log.add('reduceDimensionality', struct('method','UMAP','level','Unit'), '');
%     log.print();
%     T = log.toTable();
%     log.save('analysis_log.json');
%
%   Singleton pattern: use AnalysisLog.instance() to get a shared log.

    properties (SetAccess = private)
        Entries     struct = struct('timestamp', {}, 'method', {}, 'parameters', {}, ...
                                   'summary', {}, 'duration_s', {}, 'git_hash', {})
        SessionID   string
        StartTime   datetime
    end

    methods

        function al = AnalysisLog()
        % Constructor. Creates a new log with a unique session ID.
            al.SessionID = string(dec2hex(randi(2^32), 8));
            al.StartTime = datetime('now');
        end

        function add(al, method_name, parameters, summary, duration_s)
        % ADD  Record an analysis step.
        %
        % INPUTS:
        %   method_name - string name of the method/function called
        %   parameters  - struct of parameters used
        %   summary     - string summary of results (optional)
        %   duration_s  - elapsed time in seconds (optional)
            arguments
                al AnalysisLog
                method_name string
                parameters struct = struct()
                summary string = ""
                duration_s (1,1) double = 0
            end

            entry.timestamp = datetime('now');
            entry.method = method_name;
            entry.parameters = parameters;
            entry.summary = summary;
            entry.duration_s = duration_s;
            entry.git_hash = AnalysisLog.getGitHash();

            al.Entries(end+1) = entry;
        end

        function print(al)
        % PRINT  Display the log to the console.
            fprintf('\n=== AnalysisLog [%s] started %s ===\n', al.SessionID, string(al.StartTime));
            for i = 1:length(al.Entries)
                e = al.Entries(i);
                time_str = string(e.timestamp, 'HH:mm:ss');
                if e.duration_s > 0
                    dur_str = sprintf(' (%.1fs)', e.duration_s);
                else
                    dur_str = '';
                end
                fprintf('  [%s]%s %s', time_str, dur_str, e.method);
                if strlength(e.summary) > 0
                    fprintf(' — %s', e.summary);
                end
                fprintf('\n');
            end
            fprintf('=== %d entries ===\n\n', length(al.Entries));
        end

        function T = toTable(al)
        % TOTABLE  Export log entries as a MATLAB table.
            if isempty(al.Entries)
                T = table();
                return
            end
            timestamps = [al.Entries.timestamp]';
            methods = [al.Entries.method]';
            summaries = [al.Entries.summary]';
            durations = [al.Entries.duration_s]';
            git_hashes = string({al.Entries.git_hash})';

            % Serialize parameters to JSON strings
            param_strs = strings(length(al.Entries), 1);
            for i = 1:length(al.Entries)
                param_strs(i) = string(jsonencode(al.Entries(i).parameters));
            end

            T = table(timestamps, methods, param_strs, summaries, durations, git_hashes, ...
                'VariableNames', {'Timestamp', 'Method', 'Parameters', 'Summary', 'Duration_s', 'GitHash'});
        end

        function save(al, filepath)
        % SAVE  Export the log to a JSON file.
            arguments
                al AnalysisLog
                filepath string
            end
            log_struct.SessionID = al.SessionID;
            log_struct.StartTime = string(al.StartTime);
            log_struct.Entries = al.Entries;
            for i = 1:length(log_struct.Entries)
                log_struct.Entries(i).timestamp = string(log_struct.Entries(i).timestamp);
            end
            json_str = jsonencode(log_struct, 'PrettyPrint', true);
            fid = fopen(filepath, 'w');
            fprintf(fid, '%s', json_str);
            fclose(fid);
            fprintf('Log saved to %s\n', filepath);
        end

    end

    methods (Static)

        function al = instance()
        % INSTANCE  Get the singleton AnalysisLog (shared across all code in the session).
            persistent shared_log
            if isempty(shared_log) || ~isvalid(shared_log)
                shared_log = AnalysisLog();
            end
            al = shared_log;
        end

        function hash = getGitHash()
        % GETGITHASH  Return the current git commit hash (short), or "" if not in a repo.
            try
                [status, result] = system('git rev-parse --short HEAD');
                if status == 0
                    hash = strtrim(string(result));
                else
                    hash = "";
                end
            catch
                hash = "";
            end
        end

    end
end

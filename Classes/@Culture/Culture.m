classdef Culture < handle
% CULTURE  A set of MEArecording objects representing the same neural network
%          recorded at different conditions (e.g. drug concentrations, time points).
%
% Constructed via RecordingGroup.groupCultures(grouping_var).
%
%   culture.Recordings      - (1 x N) MEArecording array
%   culture.GroupingVar     - metadata field that differs across recordings (e.g. "Concentration")
%   culture.GroupingValues  - values of GroupingVar for each recording
%   culture.Units           - Unit array from Recordings(1) — baseline
%   culture.Metadata        - metadata struct from Recordings(1)
%   culture.N               - number of recordings in this culture

    properties
        Recordings      % (1 x N) MEArecording array
        GroupingVar     % string — metadata field that varies across recordings
        GroupingValues  % (1 x N) values of GroupingVar for each recording
    end

    properties (Dependent)
        Units           % Unit array from Recordings(1) — baseline recording
        Metadata        % metadata struct from Recordings(1)
        N               % number of recordings in this culture
    end

    methods
        function culture = Culture(recordings, grouping_var, grouping_values)
            arguments
                recordings    MEArecording = MEArecording.empty
                grouping_var  string       = ""
                grouping_values            = []
            end
            if isempty(recordings)
                return
            end
            culture.Recordings     = recordings;
            culture.GroupingVar    = grouping_var;
            culture.GroupingValues = grouping_values;
        end

        function units = get.Units(culture)
            units = culture.Recordings(1).Units;
        end

        function metadata = get.Metadata(culture)
            metadata = culture.Recordings(1).Metadata;
        end

        function n = get.N(culture)
            n = length(culture.Recordings);
        end
    end
end

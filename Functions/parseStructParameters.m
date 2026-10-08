function merged = parseStructParameters(defaults, overrides)
% PARSESTRUCTPARAMETERS  Merge an overrides struct into a defaults struct.
%
% INPUTS:
%   defaults  - Struct of default parameter values (nested structs supported)
%   overrides - Struct of user-supplied parameter overrides
%
% OUTPUT:
%   merged    - Defaults struct with overridden values applied
%
% Nested structs are merged recursively at any depth, so overriding one field of
% QC.KSLabel keeps the other QC.KSLabel defaults. A default that is an empty struct
% (e.g. Bombcell.Overrides = struct()) is a free-form slot and takes the override as is.
% Unknown fields in overrides emit a warning (parseStructParameters:unknownField) with
% their full path and are ignored. An empty overrides struct is silently ignored.

merged = mergeLevel(defaults, overrides, "");

end

function merged = mergeLevel(defaults, overrides, prefix)
    merged = defaults;
    fields = fieldnames(overrides);
    for k = 1:numel(fields)
        f = fields{k};
        path = prefix + f;
        if ~isfield(merged, f)
            warning('parseStructParameters:unknownField', 'Unrecognized parameter field: %s', path);
        elseif isstruct(overrides.(f)) && isscalar(overrides.(f)) && ...
                isstruct(merged.(f)) && isscalar(merged.(f)) && ~isempty(fieldnames(merged.(f)))
            merged.(f) = mergeLevel(merged.(f), overrides.(f), path + ".");
        else
            merged.(f) = overrides.(f);
        end
    end
end

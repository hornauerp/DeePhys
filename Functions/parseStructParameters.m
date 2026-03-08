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
% Unknown top-level or sub-level fields in overrides emit a warning.
% Empty overrides struct is silently ignored.

merged = defaults;
parameter_fields = fieldnames(overrides);

if isempty(parameter_fields)
    return
end

for pf = 1:length(parameter_fields)
    field = parameter_fields{pf};
    if ismember(field, fieldnames(merged))
        subfields = fieldnames(overrides.(field));
        for sf = 1:length(subfields)
            subfield = subfields{sf};
            if ismember(subfield, fieldnames(merged.(field)))
                merged.(field).(subfield) = overrides.(field).(subfield);
            else
                warning("Unrecognized parameter field: %s.%s", field, subfield)
            end
        end
    else
        warning("Unrecognized parameter field: %s", field)
    end
end

end

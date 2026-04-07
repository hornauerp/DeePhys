function saveDiagnosticFigure(ctc, fig, stage_name)
% SAVEDIAGNOSTICFIGURE  Save and/or hide a diagnostic figure per Diagnostics parameters.
%
%   If Diagnostics.SaveDir is non-empty, saves figure as PNG:
%     [SaveDir]/diag_[stage_name]_[yyyymmdd_HHMMSS].png
%   If Diagnostics.ShowFigures is false, closes the figure after saving.

    d = ctc.Parameters.Diagnostics;
    if ~isempty(d.SaveDir)
        if ~isfolder(d.SaveDir)
            mkdir(d.SaveDir);
        end
        ts       = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
        fname    = fullfile(d.SaveDir, sprintf('diag_%s_%s.png', stage_name, ts));
        exportgraphics(fig, fname, 'Resolution', 150);
    end
    if ~d.ShowFigures
        close(fig);
    end
end

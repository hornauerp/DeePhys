function diagnosticOptimization(ctc, results)
% DIAGNOSTICOPTIMIZATION  Diagnostic figure: "Did the optimizer converge meaningfully?"
%
% Layout (4x4 tiled figure):
%   Rows 1-2, cols 1-2: (a) Convergence trace — objective per iteration + running minimum
%   Rows 1-2, cols 3-4: (c) Best vs default   — text table with % change annotations
%   Rows 3-4, tiles 9-15: (b) Parameter marginals — one scatter per variable
%
% Called at the end of optimizeUnsupervisedUMAP when Diagnostics.Enable = true.
% Can also be called independently with its results struct.

arguments
    ctc     CellTypeClassifier
    results struct
end

% Optimised parameters: results.bestParams holds one struct per parameter group
% (UMAP, Community); defaults come from the same groups of returnDefaultParams.
def_p  = CellTypeClassifier.returnDefaultParams();
groups = fieldnames(results.bestParams);
param_names = {}; best_vals = []; def_vals = [];
for g = 1:numel(groups)
    s = results.bestParams.(groups{g});
    f = fieldnames(s);
    for k = 1:numel(f)
        param_names{end+1} = f{k}; %#ok<AGROW>
        best_vals(end+1, 1) = s.(f{k}); %#ok<AGROW>
        if isfield(def_p, groups{g}) && isfield(def_p.(groups{g}), f{k})
            def_vals(end+1, 1) = def_p.(groups{g}).(f{k}); %#ok<AGROW>
        else
            def_vals(end+1, 1) = NaN; %#ok<AGROW>
        end
    end
end
param_labels = param_names;
range_fields = strcat(param_names, 'Range');
is_integer   = ismember(param_names, {'NNeighbors', 'NDims'});
n_params     = min(numel(param_names), 7);   % marginal tiles 9..15

% Unpack BayOpt results
bo          = results.bayesoptResult;
obj_trace   = bo.ObjectiveTrace(:);          % (n_evals x 1) objective per iteration
x_trace     = bo.XTrace;                     % table: one row per eval, one col per param
n_evals     = numel(obj_trace);
iter_idx    = (1:n_evals)';

% Running minimum (ignoring NaN / failed evals)
run_min = nan(n_evals, 1);
cur_min = Inf;
for ei = 1:n_evals
    if isfinite(obj_trace(ei)) && obj_trace(ei) < cur_min
        cur_min = obj_trace(ei);
    end
    if isfinite(cur_min)
        run_min(ei) = cur_min;
    end
end

p_bo = ctc.Parameters.BayesianOptimization;

% ── Figure -------------------------------------------------------------------
fig = figure('Visible', 'on');
set(fig, 'Position', [50 50 1600 900]);
tl = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'BayOpt hyperparameter optimization diagnostics', 'FontWeight', 'bold');

% ── Tile (a): Convergence trace — rows 1-2, cols 1-2 ─────────────────────────
nexttile(1, [2, 2]);
try
    hold on;
    scatter(iter_idx, obj_trace, 30, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha', 0.5, 'DisplayName', 'Iteration');
    plot(iter_idx(isfinite(run_min)), run_min(isfinite(run_min)), ...
        '-b', 'LineWidth', 2.5, 'DisplayName', 'Running min');
    hold off;
    xlabel('Iteration');
    ylabel('Objective');
    title(sprintf('Convergence (%d evals)', n_evals), 'Interpreter', 'none');
    legend('Location', 'northeast', 'Box', 'off');
    box off;
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tile (c): Best vs default comparison — rows 1-2, cols 3-4 ────────────────
nexttile(3, [2, 2]);
try
    axis off;
    header = sprintf('%-22s  %-10s  %-10s  %-8s', 'Parameter', 'Default', 'Optimized', 'Change');
    sep    = repmat('-', 1, 56);
    lines  = {header, sep};

    for pi = 1:numel(param_names)
        def_v = def_vals(pi);
        opt_v = best_vals(pi);
        if abs(def_v) > 1e-10
            pct_change = 100 * (opt_v - def_v) / abs(def_v);
            pct_str    = sprintf('%+.0f%%', pct_change);
            if abs(pct_change) > 50
                pct_str = [pct_str ' *']; %#ok<AGROW>
            end
        else
            pct_str = 'N/A';
        end

        if is_integer(pi)
            row = sprintf('%-22s  %-10d  %-10d  %s', param_names{pi}, ...
                int32(def_v), int32(opt_v), pct_str);
        else
            row = sprintf('%-22s  %-10.4g  %-10.4g  %s', param_names{pi}, ...
                def_v, opt_v, pct_str);
        end
        lines{end+1} = row; %#ok<AGROW>
    end

    lines{end+1} = sep;
    lines{end+1} = sprintf('Best objective:         %.4f', results.bestObjective);
    lines{end+1} = '(* = >50%% change from default)';

    text(0.02, 0.97, strjoin(lines, '\n'), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontName', 'Monospaced', 'FontSize', 9, 'Interpreter', 'none');
    title('Best vs default parameters');
catch ME
    cla; axis off;
    text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 8, 'Color', [0.6 0 0]);
end

% ── Tiles (b): Parameter marginals — rows 3-4, tiles 9-15 ────────────────────
[~, best_iter] = min(obj_trace);

for pi = 1:n_params
    nexttile(8 + pi);   % tiles 9 .. 15
    try
        pname = param_names{pi};
        if ~ismember(pname, x_trace.Properties.VariableNames)
            error('"%s" not in XTrace', pname);
        end
        p_vals = x_trace.(pname);

        % Restrict to finite pairs
        valid = isfinite(obj_trace) & isfinite(p_vals);
        pv = p_vals(valid);
        ov = obj_trace(valid);

        scatter(pv, ov, 18, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.4);
        hold on;

        % Highlight best evaluation
        if best_iter <= numel(p_vals) && isfinite(p_vals(best_iter))
            scatter(p_vals(best_iter), obj_trace(best_iter), 80, 'r', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        end

        % Search-space boundaries
        if isfield(p_bo, range_fields{pi})
            range_v = p_bo.(range_fields{pi});
            xline(range_v(1), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
            xline(range_v(2), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        end
        hold off;

        xlabel(param_labels{pi}, 'FontSize', 8);
        ylabel('Obj', 'FontSize', 7);
        title(param_labels{pi}, 'FontSize', 8);
        box off;
    catch ME
        cla; axis off;
        text(0.5, 0.5, sprintf('%s\n%s', param_names{pi}, ME.message), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
end
% Tile 16 is intentionally left empty

% ── Save ─────────────────────────────────────────────────────────────────────
ctc.saveDiagnosticFigure(fig, 'optimization');

end

function plotFinalTopology(model, x, resultRoot, postResult, analyses, opts)
% PLOTFINALTOPOLOGY  Save density field and threshold topology images.
%
%   plotFinalTopology(model, x, resultRoot)
%   plotFinalTopology(model, x, resultRoot, postResult)
%   plotFinalTopology(model, x, resultRoot, postResult, analyses)
%   plotFinalTopology(..., opts)
%
%   When postResult (from postprocessSIMPResult) is supplied, an additional
%   figure comparing the best extracted topology against the raw rho>0.5
%   cut is saved as 'final_best_topology.png'.
%
%   By default this function keeps topology plotting light: it saves only the
%   rho>0.5 threshold and does not run FE metrics for threshold previews.
%   Set opts.thresholds and opts.computeMetrics to enable heavier output.

    if nargin < 4, postResult = []; end
    if nargin < 5, analyses = []; end
    if nargin < 6, opts = struct(); end
    thresholds = optField(opts, 'thresholds', 0.5);
    computeMetrics = optField(opts, 'computeMetrics', false) && ~isempty(analyses);
    saveFigFiles = optField(opts, 'saveFigFiles', false);

    fig = figure('Name', 'Final density', 'Visible', 'off');
    plotElementDensityField(model, x);
    title(sprintf('Final density, vf=%.3f', mean(x)));
    saveas(fig, fullfile(resultRoot, 'final_density.png'));
    if saveFigFiles
        savefig(fig, fullfile(resultRoot, 'final_density.fig'));
    end
    close(fig);

    xVoid = 1e-6;
    for i = 1:numel(thresholds)
        t    = thresholds(i);
        stem = sprintf('final_threshold_rho_gt_%02d', round(10 * t));
        solid_t = x > t;

        if computeMetrics
            x_bin_t = double(solid_t) + xVoid * double(~solid_t);
            perf    = evaluateStructuralPerformance(analyses, x_bin_t, 1, false);
            titleStr = sprintf('Final topology  \\rho > %.1f  |  V=%.3f  sHM=%.3e  u=%.3e', ...
                t, mean(solid_t), perf.sHM_max, perf.u_max);
        else
            titleStr = sprintf('Final topology  \\rho > %.1f', t);
        end

        fig  = figure('Name', sprintf('rho > %.1f', t), 'Visible', 'off');
        hold on; axis off; daspect([1 1 1]); view(45, 35);
        model.fe.plotSolidSelected(model.mesh.nodes, solid_t, [0.55 0.55 0.55]);
        title(titleStr);
        saveas(fig, fullfile(resultRoot, [stem '.png']));
        if saveFigFiles
            savefig(fig, fullfile(resultRoot, [stem '.fig']));
        end
        close(fig);
    end

    % Optional: comparison panel when post-processing result is available
    if isempty(postResult)
        return;
    end

    fig = figure('Name', 'Best topology comparison', 'Visible', 'off');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    hold on; axis off; daspect([1 1 1]); view(45, 35);
    model.fe.plotSolidSelected(model.mesh.nodes, x > 0.5, [0.65 0.65 0.65]);
    title(sprintf('\\rho > 0.5  (V=%.3f)', mean(x > 0.5)));

    nexttile;
    best = postResult.best;
    hold on; axis off; daspect([1 1 1]); view(45, 35);
    model.fe.plotSolidSelected(model.mesh.nodes, best.solid, [0.45 0.60 0.80]);
    if isfield(best, 'J_simp')
        metricStr = sprintf('J=%.4f', best.J_simp);
    elseif isfield(best, 'maxConstraintContinuous')
        metricStr = sprintf('maxG=%.3f', best.maxConstraintContinuous);
    elseif isfield(best, 'J_direct')
        metricStr = sprintf('J=%.4f', best.J_direct);
    else
        metricStr = '';
    end
    if isfield(best, 'sHM_max')
        title(sprintf('%s  (V=%.3f, %s, sHM=%.3e, u=%.3e)', ...
            strrep(best.label,'_',' '), best.volFrac, metricStr, best.sHM_max, best.u_max), ...
            'Interpreter', 'tex');
    else
        title(sprintf('%s  (V=%.3f, %s)', strrep(best.label,'_',' '), best.volFrac, metricStr), ...
            'Interpreter', 'tex');
    end

    sgtitle('Topology extraction: naive cut vs best method');
    saveas(fig, fullfile(resultRoot, 'final_best_topology.png'));
    if saveFigFiles
        savefig(fig, fullfile(resultRoot, 'final_best_topology.fig'));
    end
    close(fig);
end

function v = optField(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

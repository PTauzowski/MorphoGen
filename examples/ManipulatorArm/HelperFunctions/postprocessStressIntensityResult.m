function postResult = postprocessStressIntensityResult(history, model, analyses, opts)
% POSTPROCESSSTRESSINTENSITYRESULT  Post-process a stress-intensity ESO result.
%
%   postResult = postprocessStressIntensityResult(history, model, analyses, opts)
%
%   Evaluates a set of candidate ESO iterations near the target volume
%   fraction, extracts a binary topology from each, and selects the best
%   by minimising weighted compliance.  Works for both linked (testC) and
%   unlinked (testD) ESO runs.
%
%   Inputs:
%     history    struct from the ESO loop.  Required fields:
%                  .x              [nDesignVars x nIter]  density history
%                  .volumeFraction [nIter x 1]            VF at each iteration
%     model      ManipulatorModel3D  (first-config model)
%     analyses   {nConfigs x 1}  FEAnalysis objects
%     opts       struct  (all fields optional unless marked required)
%
%   Required opts fields:
%     VolFrac      target volume fraction (stopping criterion of the ESO run)
%
%   Optional opts fields:
%     expand       function handle z -> x_arm [nElems x 1]
%                  For linked runs (testC): @(z) model.segmentToArm(z)
%                  For unlinked runs (testD): @(x) x  (default = identity)
%     window       search window around VolFrac            (default: 0.05)
%     penal        SIMP exponent for compliance evaluation  (default: 3)
%     pAgg         p-norm aggregation exponent              (default: 4)
%     weights      [nConfigs x 1] per-config weights        (default: uniform)
%     C0           [nConfigs x 1] compliance normalisation
%                  (default: evaluated at the recommended best iteration)
%     runReanalysis  logical — strict binary reanalysis       (default: false)
%     useParallel  logical                                   (default: false)
%     resultRoot   path for saving outputs                   (default: '' = no saving)
%     configNames  {nConfigs x 1} string labels
%
%   Output: postResult struct with fields
%     candidates    [nCands x 1] struct  all evaluated iterations
%     best          struct  copy of the best candidate
%     bestIter      scalar  ESO iteration index of the best
%     candIters     [nCands x 1]  all evaluated iteration indices
%
%   Each candidate struct has:
%     .iter         scalar  ESO iteration index
%     .solid        [nElems x 1] logical  binary mask after cleanup
%     .volFrac      scalar
%     .J_direct     scalar  normalised p-norm compliance (near-binary field)
%     .C_direct     [nConfigs x 1]
%     .J_binary     scalar  (only if runReanalysis=true, else NaN)
%     .C_binary     [nConfigs x 1]
%     .label        string

    assert(isfield(opts, 'VolFrac'), ...
        'postprocessStressIntensityResult: opts.VolFrac is required.');

    nConfigs = numel(analyses);
    VolFrac  = opts.VolFrac;

    expand       = optField(opts, 'expand',       @(x) x);
    window       = optField(opts, 'window',       0.05);
    penal        = optField(opts, 'penal',        3);
    pAgg         = optField(opts, 'pAgg',         4);
    weights      = optField(opts, 'weights',      ones(nConfigs, 1) / nConfigs);
    runReanalysis = optField(opts, 'runReanalysis', false);
    useParallel  = optField(opts, 'useParallel',  false);
    resultRoot   = optField(opts, 'resultRoot',   '');
    configNames  = optField(opts, 'configNames',  ...
        arrayfun(@(k) sprintf('cfg%d',k), (1:nConfigs)', 'UniformOutput', false));
    plotModels   = optField(opts, 'plotModels',   {});
    plotConfigs  = optField(opts, 'plotConfigs',  []);

    weights = weights(:);

    %% -- Anchor elements for connectivity cleanup ----------------------------
    anchorElems = findAnchorElements(analyses, model.mesh);

    %% -- Identify candidate iterations --------------------------------------
    volHistory = history.volumeFraction(:);
    [candIters, recommendedIter] = selectOptimalIteration(volHistory, VolFrac, window);

    fprintf('\n[postprocess ESO] %d candidate iterations in VF window [%.3f, %.3f]: %s\n', ...
        numel(candIters), VolFrac - window, VolFrac + window, mat2str(candIters'));
    fprintf('[postprocess ESO] Recommended iteration from VF crossing: %d (VF=%.4f)\n', ...
        recommendedIter, volHistory(recommendedIter));

    %% -- Compute C0 at recommended iteration if not provided ----------------
    x_rec = expand(history.x(:, recommendedIter));
    C0 = optField(opts, 'C0', []);
    if isempty(C0)
        fprintf('[postprocess ESO] Computing C0 at recommended iteration %d...\n', recommendedIter);
        [~, ~, ~, C0] = evaluateObjectiveAndGradient( ...
            analyses, x_rec, penal, pAgg, weights, ones(nConfigs, 1), useParallel);
        C0 = max(C0, eps);
        fprintf('[postprocess ESO] C0 = [%s]\n', num2str(C0', '%.4e '));
    end

    %% -- Build candidate struct for each selected iteration -----------------
    nCands = numel(candIters);
    candidates = cell(nCands, 1);

    for i = 1:nCands
        iter   = candIters(i);
        x_iter = expand(history.x(:, iter));

        solid  = x_iter >= 0.5;
        solid  = removeDisconnectedComponents(solid, model.mesh, anchorElems);

        label  = sprintf('iter%d_vf%.3f', iter, volHistory(iter));
        cand   = buildCandidateESO(solid, x_iter, analyses, penal, pAgg, weights, C0, ...
            runReanalysis, useParallel, label);
        cand.iter = iter;
        candidates{i} = cand;
        fprintf('[postprocess ESO] iter=%3d  VF=%.4f  J_direct=%.6f', ...
            iter, cand.volFrac, cand.J_direct);
        if runReanalysis
            fprintf('  J_binary=%.6f', cand.J_binary);
        end
        fprintf('\n');
    end

    %% -- Select best candidate ----------------------------------------------
    if runReanalysis
        Jvals = cellfun(@(c) c.J_binary, candidates);
    else
        Jvals = cellfun(@(c) c.J_direct, candidates);
    end
    [~, bestIdx] = min(Jvals);
    bestCandidate = candidates{bestIdx};
    fprintf('[postprocess ESO] Best: %s  (J=%.6f, V=%.4f)\n', ...
        bestCandidate.label, Jvals(bestIdx), bestCandidate.volFrac);

    %% -- Assemble output ----------------------------------------------------
    postResult.candidates  = candidates;
    postResult.best        = bestCandidate;
    postResult.bestIter    = bestCandidate.iter;
    postResult.candIters   = candIters;

    %% -- Save to disk -------------------------------------------------------
    if ~isempty(resultRoot)
        unwrappedMode = "allSegments";
        if size(history.x, 1) == model.halfSegmentNelems
            unwrappedMode = "linked";
        end
        nElemsArm = analyses{1}.getTotalElemsNumber();
        fullPipeMetrics = evaluateStructuralPerformance(analyses, ones(nElemsArm, 1), 1, useParallel);
        saveESOPostprocessResults(model, candidates, resultRoot, plotModels, plotConfigs, ...
            unwrappedMode, fullPipeMetrics);
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function cand = buildCandidateESO(solid, x_nearBinary, analyses, penal, pAgg, ...
        weights, C0, runReanalysis, useParallel, label)

    nConfigs = numel(analyses);
    cand.label   = label;
    cand.solid   = solid;
    cand.volFrac = mean(solid);

    % Compliance at the near-binary field (fast; effectively binary for ESO)
    [J_d, ~, ~, C_d] = evaluateObjectiveAndGradient( ...
        analyses, x_nearBinary, penal, pAgg, weights, C0, useParallel);
    cand.J_direct = J_d;
    cand.C_direct = C_d;

    % Structural performance metrics on the binary topology
    xVoid_perf   = 1e-6;
    x_bin_perf   = double(solid) + (~solid) * xVoid_perf;
    perfMetrics  = evaluateStructuralPerformance(analyses, x_bin_perf, 1, useParallel);
    cand.sHM_max = perfMetrics.sHM_max;
    cand.u_max   = perfMetrics.u_max;
    cand.uz_max  = perfMetrics.uz_max;

    % Binary reanalysis (optional)
    cand.J_binary = NaN;
    cand.C_binary = NaN(nConfigs, 1);
    if runReanalysis
        xVoid  = 1e-6;
        x_bin  = double(solid) + (~solid) * xVoid;
        [J_b, ~, ~, C_b] = evaluateObjectiveAndGradient( ...
            analyses, x_bin, 1, pAgg, weights, C0, useParallel);
        cand.J_binary = J_b;
        cand.C_binary = C_b;
    end
end

% -------------------------------------------------------------------------
function saveESOPostprocessResults(model, candidates, resultRoot, plotModels, plotConfigs, ...
        unwrappedMode, fullPipeMetrics)

    % Summary CSV (full-pipe reference as first row)
    fpRow.label       = "full_pipe";
    fpRow.iter        = NaN;
    fpRow.volFrac     = 1.0;
    fpRow.J_direct    = NaN;
    fpRow.J_binary    = NaN;
    fpRow.sHM_max     = fullPipeMetrics.sHM_max;
    fpRow.u_max       = fullPipeMetrics.u_max;
    fpRow.uz_max      = fullPipeMetrics.uz_max;
    fpRow.nSolidElems = NaN;

    nCands = numel(candidates);
    rows   = cell(nCands + 1, 1);
    rows{1} = fpRow;
    for i = 1:nCands
        c = candidates{i};
        row.label       = string(c.label);
        row.iter        = c.iter;
        row.volFrac     = c.volFrac;
        row.J_direct    = c.J_direct;
        row.J_binary    = c.J_binary;
        row.sHM_max     = c.sHM_max;
        row.u_max       = c.u_max;
        row.uz_max      = c.uz_max;
        row.nSolidElems = sum(c.solid);
        rows{i + 1} = row;
    end
    T = struct2table(vertcat(rows{:}));
    writetable(T, fullfile(resultRoot, 'postprocess_eso_summary.csv'));
    fprintf('[postprocess ESO] Saved postprocess_eso_summary.csv\n');

    % Topology plot per candidate
    for i = 1:nCands
        c    = candidates{i};
        stem = sprintf('postprocess_eso_%s', c.label);
        fig  = figure('Visible', 'off', 'Name', c.label);
        hold on; axis off; daspect([1 1 1]); view(45, 35);
        model.fe.plotSolidSelected(model.mesh.nodes, c.solid, [0.45 0.60 0.80]);
        title(sprintf('%s  |  V=%.3f  J=%.4f  sHM=%.3e  uz=%.3e', ...
            strrep(c.label,'_',' '), c.volFrac, c.J_direct, c.sHM_max, c.uz_max), ...
            'Interpreter', 'tex');
        saveas(fig, fullfile(resultRoot, [stem '.png']));
        savefig(fig, fullfile(resultRoot, [stem '.fig']));
        close(fig);

        saveArmUnwrappedTopology(model, double(c.solid), resultRoot, ...
            [stem '_unwrapped'], ...
            sprintf('%s  |  V=%.3f', strrep(c.label,'_',' '), c.volFrac), ...
            struct('threshold', 0.5, 'mode', unwrappedMode));

        if ~isempty(plotModels)
            plotTopologyConfigurations(plotModels, c.solid, resultRoot, ...
                string(stem) + "_by_config", ...
                sprintf('%s by configuration', strrep(c.label, '_', ' ')), plotConfigs);
        end
    end

    Jvals = cellfun(@(c) c.J_direct, candidates);
    [~, bestIdx] = min(Jvals);
    writeTopologyPerformanceComparisonCsv(resultRoot, ...
        'postprocess_eso_topology_vs_full_pipe.csv', fullPipeMetrics, candidates, ...
        string(candidates{bestIdx}.label));

    % J vs VF comparison bar chart (if more than one candidate)
    if nCands > 1
        labels   = cellfun(@(c) c.label,   candidates, 'UniformOutput', false);
        volFracs = cellfun(@(c) c.volFrac, candidates);
        Jvals    = cellfun(@(c) c.J_direct, candidates);

        fig = figure('Visible', 'off', 'Name', 'ESO candidates');
        yyaxis left;
        bar(Jvals, 'FaceColor', [0.45 0.60 0.80]);
        ylabel('J_{direct}');
        yyaxis right;
        plot(1:nCands, volFracs, 'o-r', 'LineWidth', 1.5);
        ylabel('Volume fraction');
        set(gca, 'XTickLabel', labels, 'XTick', 1:nCands);
        xtickangle(30);
        title('ESO post-processing: candidates');
        saveas(fig, fullfile(resultRoot, 'postprocess_eso_candidates.png'));
        savefig(fig, fullfile(resultRoot, 'postprocess_eso_candidates.fig'));
        close(fig);
    end
end

% -------------------------------------------------------------------------
function v = optField(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

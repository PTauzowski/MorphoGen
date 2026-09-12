% stage2D_nonShellTopology.m
% Stage 2D: identify conditions for non-shell topology.
%
% Question: is torsion (Ms) the dominant factor enforcing closed-shell topology?
%
% Three cases:
%   A) My + Mz only (no torsion),   vf=0.40, Rfilter ~= 3 FE sizes
%   B) My only,                      vf=0.30, Rfilter ~= 3 FE sizes
%   C) My + Mz + Ms (with torsion),  vf=0.30, Rfilter ~= 3 FE sizes
%
% For each case the script saves rho, thresholded topology (rho>0.6), and
% reports: freeRhoStd, nConnComp (connected components at rho>0.6), rhoGt08.

clear; close all; clc;
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(42, 'twister');

% ----- Model geometry (from coupledExamples.m) ---------------------------
arm = armModelDefaults("thin");
E = arm.E;
nu = arm.nu;
R = arm.R;
r = arm.r;
h = arm.h;
alpha_deg = arm.alpha_deg;
res_th = arm.res_th;
wallThickness = R - r;
Rfilter = arm.Rfilter;

% ----- Optimization controls ----------------------------------------------
penal              = 4.0;
pAgg               = 4.0;   % p-norm aggregation for multi-load cases
loadTol            = 0.02;
maxIterations      = 140;
minIterations      = 25;
changeTol          = arm.mma.changeTol;
objectiveTol       = 5.0e-3;
volumeTol          = 5.0e-3;
moveLimit          = arm.mma.moveLimit;
minMoveLimit       = arm.mma.minMoveLimit;
moveDecay          = arm.mma.moveDecay;
mmaDamping         = arm.mma.mmaDamping;
useCompNorm        = true;
connectThreshold   = 0.6;   % threshold for topology connectivity analysis
denseThreshold     = 0.8;   % threshold for "dense material" fraction
useParallel = license('test', 'Distrib_Computing_Toolbox');

% ----- Case definitions ---------------------------------------------------
%   mode:        'single' | 'multi'
%   loadNames:   string array
%   loadCases:   cell array of structs
%   volFrac:     volume fraction target
%   Rfilter:     density filter radius, fixed to about 3 nominal FE sizes
cases = {
    struct('label', 'A_My_Mz', ...
           'mode',  'multi', ...
           'loadNames', {{'My','Mz'}}, ...
           'loadCases', {{struct('My',1.0), struct('Mz',1.0)}}, ...
           'volFrac', 0.40, 'Rfilter', Rfilter);
    struct('label', 'B_My', ...
           'mode',  'single', ...
           'loadNames', {{'My'}}, ...
           'loadCases', {{struct('My',1.0)}}, ...
           'volFrac', 0.30, 'Rfilter', Rfilter);
    struct('label', 'C_My_Mz_Ms', ...
           'mode',  'multi', ...
           'loadNames', {{'My','Mz','Ms'}}, ...
           'loadCases', {{struct('My',1.0), struct('Mz',1.0), struct('Ms',1.0)}}, ...
           'volFrac', 0.30, 'Rfilter', Rfilter);
};

resultRoot = fullfile(scriptDir, 'results', 'stage2D_nonShellTopology');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

fprintf('Stage 2D – non-shell topology identification\n');
fprintf('Results root: %s\n\n', resultRoot);
fprintf('Parallel outer case sweep: %d\n\n', useParallel);

summaryRows = cell(numel(cases), 1);

if useParallel
    parfor ci = 1:numel(cases)
        cs = cases{ci};
        summaryRows{ci} = runStage2DCase(cs, ci, numel(cases), resultRoot, wallThickness, ...
            E, nu, r, R, h, alpha_deg, res_th, loadTol, penal, pAgg, useCompNorm, ...
            maxIterations, minIterations, changeTol, objectiveTol, volumeTol, ...
            moveLimit, minMoveLimit, moveDecay, mmaDamping, connectThreshold, denseThreshold);
    end
else
    for ci = 1:numel(cases)
        cs = cases{ci};
        summaryRows{ci} = runStage2DCase(cs, ci, numel(cases), resultRoot, wallThickness, ...
            E, nu, r, R, h, alpha_deg, res_th, loadTol, penal, pAgg, useCompNorm, ...
            maxIterations, minIterations, changeTol, objectiveTol, volumeTol, ...
            moveLimit, minMoveLimit, moveDecay, mmaDamping, connectThreshold, denseThreshold);
    end
end

% ----- Summary table ------------------------------------------------------
T = struct2table(vertcat(summaryRows{:}));
summaryFile = fullfile(resultRoot, 'summary.csv');
writetable(T, summaryFile);
fprintf('\nSummary:\n');
disp(T(:, {'caseName','finalVF','freeRhoStd','nConnComp','rhoGt08','converged'}));
fprintf('Saved to %s\n', summaryFile);

% ----- Comparison figure across all 3 cases --------------------------------
allResults = {};
for ci = 1:numel(summaryRows)
    r_ci = summaryRows{ci};
    matFile = fullfile(resultRoot, r_ci.caseName, 'result.mat');
    if isfile(matFile)
        allResults{end+1} = load(matFile, 'rho_opt', 'const_elems'); %#ok<AGROW>
        allResults{end}.caseName = r_ci.caseName;
    end
end

if numel(allResults) == 3
    fig = figure('Visible','off','Name','Stage2D comparison');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for ci = 1:3
        nexttile;
        hold on; axis on; daspect([1 1 1]); view(45, 35);
        rho = allResults{ci}.rho_opt;
        ce  = allResults{ci}.const_elems;
        mdl_tmp = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th);
        active = rho > connectThreshold;
        mdl_tmp.fe.plotSolidSelected(mdl_tmp.mesh.nodes, active, [0.15 0.15 0.15]);
        mdl_tmp.fe.plotSolidSelected(mdl_tmp.mesh.nodes, ce,    [0.70 0.70 0.70]);
        xlabel('x'); ylabel('y'); zlabel('z');
        cs = cases{ci};
        row = summaryRows{ci};
        title(sprintf('%s\nnComp=%d  std=%.3f  >0.8: %.1f%%', ...
            strrep(cs.label,'_','\_'), row.nConnComp, row.freeRhoStd, 100*row.rhoGt08), ...
            'FontSize', 8);
    end
    sgtitle(sprintf('Stage 2D: topology at \\rho > %.1f  (A=My+Mz, B=My, C=My+Mz+Ms)', ...
        connectThreshold));
    saveas(fig, fullfile(resultRoot, 'comparison_thresholded.png'));
    close(fig);
    fprintf('Saved comparison figure.\n');
end

% ==========================================================================
function row = runStage2DCase(cs, ci, nCases, resultRoot, wallThickness, ...
        E, nu, r, R, h, alpha_deg, res_th, loadTol, penal, pAgg, useCompNorm, ...
        maxIterations, minIterations, changeTol, objectiveTol, volumeTol, ...
        moveLimit, minMoveLimit, moveDecay, mmaDamping, connectThreshold, denseThreshold)
    cs.RminFactor = cs.Rfilter / wallThickness;
    cs.loadNames = string(cs.loadNames);
    caseName = sprintf('%s_vf%03d_rmin%03d', cs.label, ...
        round(100*cs.volFrac), round(100*cs.RminFactor));
    fprintf('[%d/%d] %s\n', ci, nCases, caseName);

    caseDir = fullfile(resultRoot, caseName);
    if ~exist(caseDir, 'dir'), mkdir(caseDir); end

    row = makeEmptyRow(caseName, cs);
    try
        mdl = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th);
        mdl.analysis.isConst = true;

        allLoadsPassed = true;
        for k = 1:numel(cs.loadCases)
            mdl.applyLoadCase(cs.loadCases{k});
            lv = mdl.validateLoadApplication(cs.loadCases{k}, loadTol);
            allLoadsPassed = allLoadsPassed && lv.passed;
        end

        const_elems  = armConstRingElementIds(mdl, arm, "reference");

        if strcmp(cs.mode, 'single')
            mdl.applyLoadCase(cs.loadCases{1});
            topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance( ...
                cs.Rfilter, mdl.analysis, penal, cs.volFrac, true);
        else
            weights = ones(numel(cs.loadCases), 1) / numel(cs.loadCases);
            topOpt = SIMP_MMA_ReferenceModuleMultiLoadCompliance( ...
                cs.Rfilter, mdl, cs.loadCases, weights, pAgg, penal, ...
                cs.volFrac, true, useCompNorm);
        end

        [topOpt, free_elems, freeInitialRho] = initDesign(topOpt, const_elems, cs.volFrac);

        if strcmp(cs.mode, 'single')
            topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
            objScale = 1.0 / max(abs(topOpt.FobjValue), eps);
            history = runSingleMMA(topOpt, free_elems, const_elems, cs.volFrac, ...
                objScale, maxIterations, minIterations, changeTol, objectiveTol, ...
                volumeTol, moveLimit, minMoveLimit, moveDecay, mmaDamping);
        else
            topOpt.initializeComplianceNormalization(topOpt.x);
            topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
            objScale = 1.0 / max(abs(topOpt.FobjValue), eps);
            history = runMultiMMA(topOpt, free_elems, const_elems, cs.volFrac, ...
                objScale, maxIterations, minIterations, changeTol, objectiveTol, ...
                volumeTol, moveLimit, minMoveLimit, moveDecay, mmaDamping);
        end

        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        topOpt.computeConstraintsAndGradient(topOpt.x);
        rho_opt = topOpt.x;

        freeRho    = rho_opt(free_elems);
        freeRhoStd = std(freeRho);

        rhoAboveConn  = rho_opt > connectThreshold;
        nConnComp     = countConnectedComponents(mdl.mesh.elems, rhoAboveConn);
        rhoGt08       = mean(rho_opt > denseThreshold);

        finalVF   = mean(rho_opt);
        finalConstr = topOpt.constrValues;
        finalChange = topOpt.change;
        converged = finalChange < changeTol || ...
            objectiveHistoryConverged(history, objectiveTol, cs.mode);

        fprintf('  vf=%.4f  freeRhoStd=%.4f  nConnComp=%d  rhoGt08=%.4f\n', ...
            finalVF, freeRhoStd, nConnComp, rhoGt08);

        saveFigures(caseDir, mdl, rho_opt, const_elems, caseName, connectThreshold);
        saveHistoryCsv(caseDir, history);

        save(fullfile(caseDir, 'result.mat'), ...
            'rho_opt', 'history', 'const_elems', 'free_elems', ...
            'freeInitialRho', 'finalVF', 'finalConstr', 'finalChange', ...
            'converged', 'freeRhoStd', 'nConnComp', 'rhoGt08', ...
            'connectThreshold', 'denseThreshold', 'allLoadsPassed', ...
            'E', 'nu', 'r', 'R', 'h', 'alpha_deg', 'res_th', 'cs');

        row.status        = "ok";
        row.converged     = converged;
        row.finalVF       = finalVF;
        row.finalConstr   = finalConstr;
        row.finalChange   = finalChange;
        row.freeRhoStd    = freeRhoStd;
        row.nConnComp     = nConnComp;
        row.rhoGt08       = rhoGt08;
        row.allLoadsPassed = allLoadsPassed;
        row.iterations    = numel(history.iteration);

    catch ME
        row.status = "failed";
        row.errorMsg = string(ME.message);
        fprintf('  FAILED: %s\n', ME.message);
        save(fullfile(caseDir, 'failed.mat'), 'ME', 'cs');
    end
end

function row = makeEmptyRow(caseName, cs)
    row.caseName     = string(caseName);
    row.mode         = string(cs.mode);
    row.loads        = strjoin(cs.loadNames, '+');
    row.volFrac      = cs.volFrac;
    row.RminFactor   = cs.RminFactor;
    row.status       = "pending";
    row.errorMsg     = "";
    row.converged    = false;
    row.finalVF      = NaN;
    row.finalConstr  = NaN;
    row.finalChange  = NaN;
    row.freeRhoStd   = NaN;
    row.nConnComp    = NaN;
    row.rhoGt08      = NaN;
    row.allLoadsPassed = false;
    row.iterations   = 0;
end

function [topOpt, free_elems, freeInitRho] = initDesign(topOpt, const_elems, volFrac)
    topOpt.setConstElems(const_elems);
    free_elems = setdiff((1:topOpt.totalFENumber)', const_elems);
    freeInitRho = (volFrac * topOpt.totalFENumber - numel(const_elems)) / numel(free_elems);
    assert(freeInitRho > topOpt.xmin(1), ...
        'Volume fraction %.3f infeasible with %d const elems.', volFrac, numel(const_elems));
    topOpt.x(:) = freeInitRho;
    topOpt.x(const_elems) = 1.0;
    topOpt.xmin(const_elems) = 1.0 - 1e-9;
    topOpt.xmax(const_elems) = 1.0;
    topOpt.xold1 = topOpt.x;
    topOpt.xold2 = topOpt.x;
end

function history = runSingleMMA(topOpt, free_elems, const_elems, volFrac, objScale, ...
        maxIter, minIter, changeTol, objTol, volTol, moveLimit, minMove, moveDecay, damp)
    history = initHistory(maxIter, 1, 'single');
    topOpt.iteration = 1;
    topOpt.resetAnalysis();
    while topOpt.iteration <= maxIter
        prev = topOpt.x;
        stepMMA(topOpt, objScale);
        topOpt.x = dampProject(topOpt.x, prev, free_elems, const_elems, volFrac, ...
            topOpt.xmin, topOpt.xmax, topOpt.iteration, moveLimit, minMove, moveDecay, damp);
        topOpt.change = max(abs(topOpt.x - prev));
        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        topOpt.computeConstraintsAndGradient(topOpt.x);
        k = topOpt.iteration;
        history.iteration(k)  = k;
        history.objective(k)  = topOpt.FobjValue;
        history.compliance(k) = topOpt.FobjValue;
        history.normalizedCompliance(k) = NaN;
        history.volumeFraction(k) = mean(topOpt.x);
        history.constraint(k) = topOpt.constrValues;
        history.change(k)     = topOpt.change;
        topOpt.iteration = topOpt.iteration + 1;
        if topOpt.iteration > minIter && ...
                (topOpt.change < changeTol || objectiveHistoryConverged(history, objTol, 'single')) && ...
                abs(topOpt.constrValues) < volTol
            break;
        end
    end
    history = trimHistory(history);
end

function history = runMultiMMA(topOpt, free_elems, const_elems, volFrac, objScale, ...
        maxIter, minIter, changeTol, objTol, volTol, moveLimit, minMove, moveDecay, damp)
    nLoads = numel(topOpt.loadCases);
    history = initHistory(maxIter, nLoads, 'multi');
    topOpt.iteration = 1;
    topOpt.resetAnalysis();
    while topOpt.iteration <= maxIter
        prev = topOpt.x;
        stepMMA(topOpt, objScale);
        topOpt.x = dampProject(topOpt.x, prev, free_elems, const_elems, volFrac, ...
            topOpt.xmin, topOpt.xmax, topOpt.iteration, moveLimit, minMove, moveDecay, damp);
        topOpt.change = max(abs(topOpt.x - prev));
        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        topOpt.computeConstraintsAndGradient(topOpt.x);
        k = topOpt.iteration;
        history.iteration(k)  = k;
        history.objective(k)  = topOpt.FobjValue;
        history.compliance(k, :) = topOpt.complianceValues(:)';
        history.normalizedCompliance(k, :) = topOpt.normalizedComplianceValues(:)';
        history.volumeFraction(k) = mean(topOpt.x);
        history.constraint(k) = topOpt.constrValues;
        history.change(k)     = topOpt.change;
        topOpt.iteration = topOpt.iteration + 1;
        if topOpt.iteration > minIter && ...
                (topOpt.change < changeTol || objectiveHistoryConverged(history, objTol, 'multi')) && ...
                abs(topOpt.constrValues) < volTol
            break;
        end
    end
    history = trimHistory(history);
end

function stepMMA(topOpt, objScale)
    topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
    topOpt.computeConstraintsAndGradient(topOpt.x);
    topOpt.gradFobjValue = topOpt.filteringByMAtrix(topOpt.gradFobjValue);
    f0val = objScale * topOpt.FobjValue;
    df0dx = objScale * topOpt.gradFobjValue;
    [xmma,~,~,~,~,~,~,~,~,topOpt.low,topOpt.upp] = mmasub2( ...
        size(topOpt.constrValues, 1), topOpt.totalFENumber, topOpt.iteration, ...
        topOpt.x, topOpt.xmin, topOpt.xmax, topOpt.xold1, topOpt.xold2, ...
        f0val, df0dx, 0*df0dx, ...
        topOpt.constrValues, topOpt.gradConstrValues, 0*topOpt.gradConstrValues, ...
        topOpt.low, topOpt.upp, topOpt.a0, topOpt.ai, topOpt.ci, topOpt.di);
    if topOpt.iteration > 1, topOpt.xold2 = topOpt.xold1; end
    topOpt.xold1 = topOpt.x;
    topOpt.x = xmma;
end

function x = dampProject(x, prev, free_elems, const_elems, volFrac, xmin, xmax, ...
        iter, moveLimit, minMove, moveDecay, damp)
    ml = max(minMove, moveLimit * moveDecay^(iter - 1));
    x(free_elems) = min(max(x(free_elems), prev(free_elems)-ml), prev(free_elems)+ml);
    x(const_elems) = 1.0;
    x = prev + damp*(x - prev);
    x(const_elems) = 1.0;
    x = enforceVF(x, free_elems, const_elems, volFrac, xmin, xmax);
end

function x = enforceVF(x, free_elems, const_elems, volFrac, xmin, xmax)
    x = min(max(x, xmin), xmax);
    x(const_elems) = 1.0;
    target = volFrac*numel(x) - sum(x(const_elems));
    target = min(max(target, sum(xmin(free_elems))), sum(xmax(free_elems)));
    lo = min(xmin(free_elems) - x(free_elems));
    hi = max(xmax(free_elems) - x(free_elems));
    for k = 1:60
        shift = 0.5*(lo+hi);
        cand = min(max(x(free_elems)+shift, xmin(free_elems)), xmax(free_elems));
        if sum(cand) < target, lo = shift; else, hi = shift; end
    end
    x(free_elems) = min(max(x(free_elems)+0.5*(lo+hi), xmin(free_elems)), xmax(free_elems));
    x(const_elems) = 1.0;
end

function history = initHistory(nIter, nLoads, mode)
    history.mode = string(mode);
    history.iteration         = nan(nIter, 1);
    history.objective         = nan(nIter, 1);
    history.compliance        = nan(nIter, nLoads);
    history.normalizedCompliance = nan(nIter, nLoads);
    history.volumeFraction    = nan(nIter, 1);
    history.constraint        = nan(nIter, 1);
    history.change            = nan(nIter, 1);
end

function history = trimHistory(history)
    keep = ~isnan(history.iteration);
    fields = fieldnames(history);
    for i = 1:numel(fields)
        if isnumeric(history.(fields{i})) && size(history.(fields{i}), 1) == numel(keep)
            history.(fields{i}) = history.(fields{i})(keep, :);
        end
    end
end

function tf = objectiveHistoryConverged(history, tol, mode)
    if strcmp(mode, 'single')
        obj = history.compliance(~isnan(history.compliance(:,1)), 1);
    else
        obj = history.objective(~isnan(history.objective));
    end
    nW = min(10, numel(obj));
    if nW < 3, tf = false; return; end
    j = obj(end-nW+1:end);
    tf = (max(j)-min(j)) / max(abs(j(end)), eps) < tol;
end

function nComp = countConnectedComponents(elems, activeElems)
    % Count connected components among active elements using face adjacency.
    % Two hex8 elements are face-adjacent if they share exactly 4 nodes.
    activeIds = find(activeElems(:));
    nActive   = numel(activeIds);
    if nActive == 0
        nComp = 0;
        return;
    end
    nNodes       = max(elems(:));
    nodesPerElem = size(elems, 2);
    activeElNodes = elems(activeIds, :);               % nActive x nodesPerElem
    localIdx = repmat((1:nActive)', 1, nodesPerElem);
    % Sparse incidence: entry (node, localElem) = 1
    N2E = sparse(activeElNodes(:), localIdx(:), 1, nNodes, nActive);
    % Shared-node count: (i,j) = number of nodes shared between active elems i and j
    shared = N2E' * N2E;   % nActive x nActive sparse
    % Face adjacency for hex8: share >= 4 nodes (degenerate quads possible → keep >=4)
    [ii, jj, vv] = find(shared);
    mask = (vv >= 4) & (ii ~= jj);
    if ~any(mask)
        nComp = nActive;   % every active element is isolated
        return;
    end
    G = graph(ii(mask), jj(mask), [], nActive);
    cc = conncomp(G);
    nComp = max(cc);
end

function saveFigures(caseDir, mdl, rho, const_elems, caseName, threshold)
    titleStr = strrep(caseName, '_', '\_');

    % Density field
    fig = figure('Visible','off','Name',[caseName ' density']);
    plotDensityField(mdl, rho, const_elems);
    title([titleStr ' – density']); colorbar; caxis([0 1]);
    saveas(fig, fullfile(caseDir, 'density_field.png'));
    savefig(fig, fullfile(caseDir, 'density_field.fig'));
    close(fig);

    % Density histogram
    fig = figure('Visible','off','Name',[caseName ' histogram']);
    histogram(rho, 20, 'BinLimits',[0 1]); grid on;
    xlabel('\rho'); ylabel('Element count');
    title([titleStr ' – histogram']);
    saveas(fig, fullfile(caseDir, 'density_histogram.png'));
    savefig(fig, fullfile(caseDir, 'density_histogram.fig'));
    close(fig);

    stem_t = sprintf('topology_rho_gt_%02d', round(10*threshold));
    plotTopology(mdl, rho, caseDir, struct('smoothed', true, 'saveFig', true, 'threshold', threshold, 'filenameStem', stem_t));
    exportTopology(mdl, rho, caseDir, struct('threshold', threshold, 'filenameStem', stem_t));

    plotTopology(mdl, rho, caseDir, struct('smoothed', true, 'saveFig', true, 'threshold', 0.8, 'filenameStem', 'topology_rho_gt_08'));
    exportTopology(mdl, rho, caseDir, struct('threshold', 0.8, 'filenameStem', 'topology_rho_gt_08'));
end

function saveHistoryCsv(caseDir, history)
    T = table(history.iteration(:), history.objective(:), ...
        history.volumeFraction(:), history.constraint(:), history.change(:), ...
        'VariableNames', {'iteration', 'objective', 'volumeFraction', 'constraint', 'change'});

    if isfield(history, 'compliance')
        T = [T array2table(history.compliance, ...
            'VariableNames', numberedNames('compliance', size(history.compliance, 2)))]; %#ok<AGROW>
    end
    if isfield(history, 'normalizedCompliance')
        T = [T array2table(history.normalizedCompliance, ...
            'VariableNames', numberedNames('normalizedCompliance', size(history.normalizedCompliance, 2)))]; %#ok<AGROW>
    end

    writetable(T, fullfile(caseDir, 'history.csv'));
end

function names = numberedNames(prefix, n)
    names = arrayfun(@(i) sprintf('%s_%d', prefix, i), 1:n, 'UniformOutput', false);
end

function plotDensityField(mdl, rho, const_elems)
    hold on; axis on; daspect([1 1 1]); view(45, 35);
    colormap(parula);
    elems = mdl.mesh.elems;
    facePattern = mdl.fe.sf.fcontours';
    nFacesPerElem = size(facePattern, 1);
    faces = zeros(size(elems,1)*nFacesPerElem, size(facePattern,2));
    faceColor = zeros(size(faces,1), 1);
    row = 1;
    for e = 1:size(elems,1)
        ef = reshape(elems(e, facePattern(:)), size(facePattern));
        n  = size(ef, 1);
        faces(row:row+n-1, :) = ef;
        faceColor(row:row+n-1) = rho(e);
        row = row + n;
    end
    patch('Vertices', mdl.mesh.nodes, 'Faces', faces, ...
        'FaceVertexCData', faceColor, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 1.0);
    mdl.fe.plotSolidSelected(mdl.mesh.nodes, const_elems, [0.70 0.70 0.70]);
    xlabel('x'); ylabel('y'); zlabel('z');
end

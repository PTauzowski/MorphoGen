% sweepReferenceModuleSIMP
% Stage 2C diagnostic sweep for ReferenceModuleSolidModel SIMP topologies.
%
% This script intentionally keeps the existing load application and
% compliance sensitivity implementations unchanged. It only varies load
% combinations, volume fraction, filter radius, and p-norm aggregation.

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(4, 'twister');

% ----- Model parameters from coupledExamples.m ---------------------------
E = 2.0e9;
nu = 0.35;
R = 0.14;
r = 0.08;
alpha = 22.5;
segmentLength = 0.25;
res_thickness = 4;

h = segmentLength;
alpha_deg = alpha;
res_th = res_thickness;
wallThickness = R - r;

% ----- Sweep parameters --------------------------------------------------
thresholdValues = [0.5, 0.7];

caseSpecs = {
    makeCaseParams('single', 'My_only', string({'My'}), ...
        {struct('My', 1.0)}, 0.4, 0.5, NaN, wallThickness);
    makeCaseParams('single', 'Ms_only', string({'Ms'}), ...
        {struct('Ms', 1.0)}, 0.4, 0.5, NaN, wallThickness);
    makeCaseParams('single', 'Ty_only', string({'Ty'}), ...
        {struct('Ty', 1.0)}, 0.4, 0.5, NaN, wallThickness);
    makeCaseParams('multi', 'My_Mz_Ms_Ty_Tz_baseline', string({'My','Mz','Ms','Ty','Tz'}), ...
        {struct('My', 1.0), struct('Mz', 1.0), struct('Ms', 1.0), struct('Ty', 1.0), struct('Tz', 1.0)}, 0.5, 1.5, 4, wallThickness);
    makeCaseParams('multi', 'My_Mz_Ms_Ty_Tz_vf040_rmin150', string({'My','Mz','Ms','Ty','Tz'}), ...
        {struct('My', 1.0), struct('Mz', 1.0), struct('Ms', 1.0), struct('Ty', 1.0), struct('Tz', 1.0)}, 0.4, 1.5, 4, wallThickness);
    makeCaseParams('multi', 'My_Mz_Ms_Ty_Tz_vf030_rmin150', string({'My','Mz','Ms','Ty','Tz'}), ...
        {struct('My', 1.0), struct('Mz', 1.0), struct('Ms', 1.0), struct('Ty', 1.0), struct('Tz', 1.0)}, 0.3, 1.5, 4, wallThickness);
    makeCaseParams('multi', 'My_Mz_Ms_Ty_Tz_vf040_rmin100', string({'My','Mz','Ms','Ty','Tz'}), ...
        {struct('My', 1.0), struct('Mz', 1.0), struct('Ms', 1.0), struct('Ty', 1.0), struct('Tz', 1.0)}, 0.4, 1.0, 4, wallThickness);
    makeCaseParams('multi', 'My_Mz_Ms_Ty_Tz_vf040_rmin050_p04', string({'My','Mz','Ms','Ty','Tz'}), ...
        {struct('My', 1.0), struct('Mz', 1.0), struct('Ms', 1.0), struct('Ty', 1.0), struct('Tz', 1.0)}, 0.4, 0.5, 4, wallThickness);
    makeCaseParams('multi', 'My_Mz_Ms_Ty_Tz_vf040_rmin050_p08', string({'My','Mz','Ms','Ty','Tz'}), ...
        {struct('My', 1.0), struct('Mz', 1.0), struct('Ms', 1.0), struct('Ty', 1.0), struct('Tz', 1.0)}, 0.4, 0.5, 8, wallThickness);
};

% ----- Optimization controls copied from Stage 2B/2C scripts -------------
loadTol = 0.02;
penal = 4.0;
useComplianceNormalization = true;
maxIterations = 140;
minIterations = 25;
changeTol = 1.0e-3;
objectiveTol = 5.0e-3;
volumeTol = 5.0e-3;
moveLimit = 0.05;
minMoveLimit = 0.003;
moveDecay = 0.95;
mmaDamping = 0.25;

resultRoot = fullfile(scriptDir, 'results', 'referenceModuleSIMP_sweep');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

summaryRows = {};
totalCases = numel(caseSpecs);

fprintf('Reference-module SIMP diagnostic sweep: %d cases\n', totalCases);
fprintf('Results root: %s\n\n', resultRoot);

% ----- Requested diagnostic cases ---------------------------------------
for caseIndex = 1:totalCases
    params = caseSpecs{caseIndex};
    params.resultRoot = resultRoot;
    fprintf('\n[%03d/%03d] %s\n', caseIndex, totalCases, params.caseName);
    summaryRows{end+1, 1} = runCase(params); %#ok<SAGROW>
end

summaryTable = struct2table(vertcat(summaryRows{:}));
summaryFile = fullfile(resultRoot, 'summary.csv');
writetable(summaryTable, summaryFile);
fprintf('\nSaved summary table to %s\n', summaryFile);

% =========================================================================
function params = makeCaseParams(mode, loadSetName, loadNames, loadCases, volFracTarget, RminFactor, pAgg, wallThickness)
    if strcmp(mode, 'single')
        caseName = sprintf('%s_vf%03d_rmin%03d', ...
            loadSetName, round(100 * volFracTarget), round(100 * RminFactor));
    else
        caseName = sprintf('%s_vf%03d_rmin%03d_p%02d', ...
            loadSetName, round(100 * volFracTarget), round(100 * RminFactor), round(pAgg));
    end

    params.mode = mode;
    params.loadSetName = loadSetName;
    params.loadNames = loadNames(:);
    params.loadCases = loadCases(:);
    params.volFracTarget = volFracTarget;
    params.RminFactor = RminFactor;
    params.Rfilter = RminFactor * wallThickness;
    params.pAgg = pAgg;
    params.caseName = caseName;
end

function row = runCase(params)
    caseDir = fullfile(params.resultRoot, params.caseName);
    if ~exist(caseDir, 'dir')
        mkdir(caseDir);
    end

    E = 2.0e9;
    nu = 0.35;
    R = 0.14;
    r = 0.08;
    h = 0.25;
    alpha_deg = 22.5;
    res_th = 4;

    loadTol = 0.02;
    penal = 4.0;
    useComplianceNormalization = true;
    maxIterations = 140;
    minIterations = 25;
    changeTol = 1.0e-3;
    objectiveTol = 5.0e-3;
    volumeTol = 5.0e-3;
    moveLimit = 0.05;
    minMoveLimit = 0.003;
    moveDecay = 0.95;
    mmaDamping = 0.25;
    thresholdValues = [0.5, 0.7];

    row = emptySummaryRow(params, caseDir);
    try
        mdl = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th);
        mdl.analysis.isConst = true;

        loadValidation = cell(numel(params.loadCases), 1);
        allLoadsPassed = true;
        for k = 1:numel(params.loadCases)
            mdl.applyLoadCase(params.loadCases{k});
            loadValidation{k} = mdl.validateLoadApplication(params.loadCases{k}, loadTol);
            allLoadsPassed = allLoadsPassed && loadValidation{k}.passed;
        end

        fixedNodeIds = find(mdl.fixedFaceSelector.select(mdl.mesh.nodes));
        const_elems = find(any(ismember(mdl.mesh.elems, [mdl.loaded_node_ids(:); fixedNodeIds(:)]), 2));
        const_elems = unique(const_elems(:));

        if strcmp(params.mode, 'single')
            mdl.applyLoadCase(params.loadCases{1});
            topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance( ...
                params.Rfilter, mdl.analysis, penal, params.volFracTarget, true);
            [topOpt, free_elems, freeInitialRho] = initializeDesign(topOpt, const_elems, params.volFracTarget);
            topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
            complianceNormalization = max(topOpt.FobjValue, eps);
            objectiveScale = 1.0 / max(abs(topOpt.FobjValue), eps);
            history = runSingleLoadMMA(topOpt, free_elems, const_elems, params.volFracTarget, ...
                objectiveScale, maxIterations, minIterations, changeTol, objectiveTol, ...
                volumeTol, moveLimit, minMoveLimit, moveDecay, mmaDamping);
            topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
            topOpt.computeConstraintsAndGradient(topOpt.x);
            C_k = topOpt.FobjValue;
            normalized_C_k = C_k / complianceNormalization;
            J = C_k;
            complianceNormalizationValues = complianceNormalization;
            weights = 1.0;
            loadContributions = 1.0;
        else
            weights = ones(numel(params.loadCases), 1) / numel(params.loadCases);
            topOpt = SIMP_MMA_ReferenceModuleMultiLoadCompliance( ...
                params.Rfilter, mdl, params.loadCases, weights, params.pAgg, penal, ...
                params.volFracTarget, true, useComplianceNormalization);
            [topOpt, free_elems, freeInitialRho] = initializeDesign(topOpt, const_elems, params.volFracTarget);
            topOpt.initializeComplianceNormalization(topOpt.x);
            topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
            objectiveScale = 1.0 / max(abs(topOpt.FobjValue), eps);
            history = runMultiLoadMMA(topOpt, free_elems, const_elems, params.volFracTarget, ...
                objectiveScale, maxIterations, minIterations, changeTol, objectiveTol, ...
                volumeTol, moveLimit, minMoveLimit, moveDecay, mmaDamping);
            topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
            topOpt.computeConstraintsAndGradient(topOpt.x);
            C_k = topOpt.complianceValues(:);
            normalized_C_k = topOpt.normalizedComplianceValues(:);
            J = topOpt.FobjValue;
            complianceNormalizationValues = topOpt.complianceNormalization(:);
            loadContributions = weights .* (normalized_C_k .^ params.pAgg);
            loadContributions = loadContributions / sum(loadContributions);
        end

        rho_opt = topOpt.x;
        finalVolumeFraction = mean(rho_opt);
        finalConstraint = topOpt.constrValues;
        finalChange = topOpt.change;
        mmaConverged = finalChange < changeTol || objectiveHistoryConverged(history, objectiveTol, params.mode);
        volumeActive = abs(finalConstraint) < volumeTol;
        freeRho = rho_opt(free_elems);
        freeRhoStd = std(freeRho);
        freeRhoRange = max(freeRho) - min(freeRho);
        densityFractions = arrayfun(@(t) mean(rho_opt > t), thresholdValues);
        densityHistogram.edges = linspace(0, 1, 21);
        densityHistogram.counts = histcounts(rho_opt, densityHistogram.edges);

        saveCaseOutputs(caseDir, mdl, rho_opt, const_elems, thresholdValues, params.caseName);

        resultFile = fullfile(caseDir, 'result.mat');
        save(resultFile, ...
            'rho_opt', 'history', 'C_k', 'normalized_C_k', 'J', ...
            'complianceNormalizationValues', 'loadContributions', 'loadValidation', ...
            'allLoadsPassed', 'const_elems', 'free_elems', 'freeInitialRho', ...
            'finalVolumeFraction', 'finalConstraint', 'finalChange', 'mmaConverged', ...
            'volumeActive', 'freeRhoStd', 'freeRhoRange', 'densityFractions', ...
            'densityHistogram', 'params', 'weights', 'penal', 'objectiveScale', 'moveLimit', ...
            'minMoveLimit', 'moveDecay', 'mmaDamping', 'E', 'nu', 'r', 'R', ...
            'h', 'alpha_deg', 'res_th');

        row.status = "ok";
        row.resultFile = string(resultFile);
        row.allLoadsPassed = allLoadsPassed;
        row.iterations = numel(history.iteration);
        row.mmaConverged = mmaConverged;
        row.volumeActive = volumeActive;
        row.finalVolumeFraction = finalVolumeFraction;
        row.finalConstraint = finalConstraint;
        row.finalChange = finalChange;
        row.J = J;
        row.Cmax = max(C_k);
        row.Cmean = mean(C_k);
        row.normCmax = max(normalized_C_k);
        row.normCmean = mean(normalized_C_k);
        row.freeRhoStd = freeRhoStd;
        row.freeRhoRange = freeRhoRange;
        row.rhoGt05 = densityFractions(1);
        row.rhoGt07 = densityFractions(2);
        row.maxLoadContribution = max(loadContributions);

        fprintf('  J=%.6e, vf=%.4f, change=%.3e, status=ok\n', J, finalVolumeFraction, finalChange);
    catch ME
        row.status = "failed";
        row.errorMessage = string(ME.message);
        save(fullfile(caseDir, 'failed.mat'), 'ME', 'params');
        fprintf('  FAILED: %s\n', ME.message);
    end
end

function row = emptySummaryRow(params, caseDir)
    row.caseName = string(params.caseName);
    row.mode = string(params.mode);
    row.loadSet = string(params.loadSetName);
    row.loads = strjoin(params.loadNames, '+');
    row.volFracTarget = params.volFracTarget;
    row.RminFactor = params.RminFactor;
    row.Rfilter = params.Rfilter;
    row.pAgg = params.pAgg;
    row.status = "pending";
    row.errorMessage = "";
    row.caseDir = string(caseDir);
    row.resultFile = "";
    row.allLoadsPassed = false;
    row.iterations = 0;
    row.mmaConverged = false;
    row.volumeActive = false;
    row.finalVolumeFraction = NaN;
    row.finalConstraint = NaN;
    row.finalChange = NaN;
    row.J = NaN;
    row.Cmax = NaN;
    row.Cmean = NaN;
    row.normCmax = NaN;
    row.normCmean = NaN;
    row.freeRhoStd = NaN;
    row.freeRhoRange = NaN;
    row.rhoGt05 = NaN;
    row.rhoGt07 = NaN;
    row.maxLoadContribution = NaN;
end

function [topOpt, free_elems, freeInitialRho] = initializeDesign(topOpt, const_elems, volFracTarget)
    topOpt.setConstElems(const_elems);
    free_elems = setdiff((1:topOpt.totalFENumber)', const_elems);
    freeInitialRho = (volFracTarget * topOpt.totalFENumber - numel(const_elems)) / numel(free_elems);
    assert(freeInitialRho > topOpt.xmin(1), ...
        'Requested volume fraction %.3f is infeasible with %d constant elements.', ...
        volFracTarget, numel(const_elems));

    topOpt.x(:) = freeInitialRho;
    topOpt.x(const_elems) = 1.0;
    topOpt.xmin(const_elems) = 1.0 - 1.0e-9;
    topOpt.xmax(const_elems) = 1.0;
    topOpt.xold1 = topOpt.x;
    topOpt.xold2 = topOpt.x;
end

function history = runSingleLoadMMA(topOpt, free_elems, const_elems, volFracTarget, objectiveScale, ...
        maxIterations, minIterations, changeTol, objectiveTol, volumeTol, ...
        moveLimit, minMoveLimit, moveDecay, mmaDamping)
    history = initializeHistory(maxIterations, 1, 'single');
    topOpt.iteration = 1;
    topOpt.resetAnalysis();

    while topOpt.iteration <= maxIterations
        previousRho = topOpt.x;
        updateDesignScaledMMA(topOpt, objectiveScale);
        topOpt.x = dampAndProject(topOpt.x, previousRho, free_elems, const_elems, ...
            volFracTarget, topOpt.xmin, topOpt.xmax, topOpt.iteration, ...
            moveLimit, minMoveLimit, moveDecay, mmaDamping);
        topOpt.change = max(abs(topOpt.x - previousRho));
        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        topOpt.computeConstraintsAndGradient(topOpt.x);

        k = topOpt.iteration;
        history.iteration(k, 1) = k;
        history.objective(k, 1) = topOpt.FobjValue;
        history.compliance(k, 1) = topOpt.FobjValue;
        history.normalizedCompliance(k, 1) = NaN;
        history.volumeFraction(k, 1) = mean(topOpt.x);
        history.constraint(k, 1) = topOpt.constrValues;
        history.change(k, 1) = topOpt.change;

        topOpt.iteration = topOpt.iteration + 1;
        if topOpt.iteration > minIterations && ...
                (topOpt.change < changeTol || objectiveHistoryConverged(history, objectiveTol, 'single')) && ...
                abs(topOpt.constrValues) < volumeTol
            break;
        end
    end
    history = trimHistory(history);
end

function history = runMultiLoadMMA(topOpt, free_elems, const_elems, volFracTarget, objectiveScale, ...
        maxIterations, minIterations, changeTol, objectiveTol, volumeTol, ...
        moveLimit, minMoveLimit, moveDecay, mmaDamping)
    history = initializeHistory(maxIterations, numel(topOpt.loadCases), 'multi');
    topOpt.iteration = 1;
    topOpt.resetAnalysis();

    while topOpt.iteration <= maxIterations
        previousRho = topOpt.x;
        updateDesignScaledMMA(topOpt, objectiveScale);
        topOpt.x = dampAndProject(topOpt.x, previousRho, free_elems, const_elems, ...
            volFracTarget, topOpt.xmin, topOpt.xmax, topOpt.iteration, ...
            moveLimit, minMoveLimit, moveDecay, mmaDamping);
        topOpt.change = max(abs(topOpt.x - previousRho));
        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        topOpt.computeConstraintsAndGradient(topOpt.x);

        k = topOpt.iteration;
        history.iteration(k, 1) = k;
        history.objective(k, 1) = topOpt.FobjValue;
        history.compliance(k, :) = topOpt.complianceValues(:)';
        history.normalizedCompliance(k, :) = topOpt.normalizedComplianceValues(:)';
        history.volumeFraction(k, 1) = mean(topOpt.x);
        history.constraint(k, 1) = topOpt.constrValues;
        history.change(k, 1) = topOpt.change;

        topOpt.iteration = topOpt.iteration + 1;
        if topOpt.iteration > minIterations && ...
                (topOpt.change < changeTol || objectiveHistoryConverged(history, objectiveTol, 'multi')) && ...
                abs(topOpt.constrValues) < volumeTol
            break;
        end
    end
    history = trimHistory(history);
end

function x = dampAndProject(candidateRho, previousRho, free_elems, const_elems, volFracTarget, ...
        xmin, xmax, iteration, moveLimit, minMoveLimit, moveDecay, mmaDamping)
    currentMoveLimit = max(minMoveLimit, moveLimit * moveDecay^(iteration - 1));
    candidateRho(free_elems) = min(max(candidateRho(free_elems), ...
        previousRho(free_elems) - currentMoveLimit), previousRho(free_elems) + currentMoveLimit);
    candidateRho(const_elems) = 1.0;
    candidateRho = previousRho + mmaDamping * (candidateRho - previousRho);
    candidateRho(const_elems) = 1.0;
    x = enforceVolumeFraction(candidateRho, free_elems, const_elems, volFracTarget, xmin, xmax);
end

function history = initializeHistory(nIter, nLoads, mode)
    history.mode = string(mode);
    history.iteration = nan(nIter, 1);
    history.objective = nan(nIter, 1);
    history.compliance = nan(nIter, nLoads);
    history.normalizedCompliance = nan(nIter, nLoads);
    history.volumeFraction = nan(nIter, 1);
    history.constraint = nan(nIter, 1);
    history.change = nan(nIter, 1);
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
        objective = history.compliance(~isnan(history.compliance(:, 1)), 1);
    else
        objective = history.objective(~isnan(history.objective));
    end
    nWindow = min(10, numel(objective));
    if nWindow < 3
        tf = false;
        return;
    end
    j = objective(end-nWindow+1:end);
    tf = (max(j) - min(j)) / max(abs(j(end)), eps) < tol;
end

function updateDesignScaledMMA(topOpt, objectiveScale)
    topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
    topOpt.computeConstraintsAndGradient(topOpt.x);
    topOpt.gradFobjValue = topOpt.filteringByMAtrix(topOpt.gradFobjValue);

    f0val = objectiveScale * topOpt.FobjValue;
    df0dx = objectiveScale * topOpt.gradFobjValue;

    [xmma,~,~,~,~,~,~,~,~,topOpt.low,topOpt.upp] = mmasub2( ...
        size(topOpt.constrValues, 1), ...
        topOpt.totalFENumber, ...
        topOpt.iteration, ...
        topOpt.x, ...
        topOpt.xmin, topOpt.xmax, ...
        topOpt.xold1, topOpt.xold2, ...
        f0val, df0dx, 0 * df0dx, ...
        topOpt.constrValues, topOpt.gradConstrValues, 0 * topOpt.gradConstrValues, ...
        topOpt.low, topOpt.upp, topOpt.a0, topOpt.ai, topOpt.ci, topOpt.di);

    if topOpt.iteration > 1
        topOpt.xold2 = topOpt.xold1;
    end
    topOpt.xold1 = topOpt.x;
    topOpt.x = xmma;
end

function x = enforceVolumeFraction(x, free_elems, const_elems, volFracTarget, xmin, xmax)
    x = min(max(x, xmin), xmax);
    x(const_elems) = 1.0;
    targetFreeVolume = volFracTarget * numel(x) - sum(x(const_elems));
    targetFreeVolume = min(max(targetFreeVolume, sum(xmin(free_elems))), sum(xmax(free_elems)));

    lo = min(xmin(free_elems) - x(free_elems));
    hi = max(xmax(free_elems) - x(free_elems));
    for k = 1:60
        shift = 0.5 * (lo + hi);
        candidate = min(max(x(free_elems) + shift, xmin(free_elems)), xmax(free_elems));
        if sum(candidate) < targetFreeVolume
            lo = shift;
        else
            hi = shift;
        end
    end

    x(free_elems) = min(max(x(free_elems) + 0.5 * (lo + hi), xmin(free_elems)), xmax(free_elems));
    x(const_elems) = 1.0;
end

function saveCaseOutputs(caseDir, mdl, rho, const_elems, thresholdValues, caseName)
    fig = figure('Visible', 'off', 'Name', [caseName ' density']);
    plotElementDensityField(mdl, rho, const_elems);
    title(sprintf('%s density', strrep(caseName, '_', '\_')));
    saveas(fig, fullfile(caseDir, 'final_topology.png'));
    close(fig);

    fig = figure('Visible', 'off', 'Name', [caseName ' histogram']);
    histogram(rho, 20, 'BinLimits', [0 1]);
    grid on; xlabel('\rho'); ylabel('Element count');
    title(sprintf('%s density histogram', strrep(caseName, '_', '\_')));
    saveas(fig, fullfile(caseDir, 'density_histogram.png'));
    close(fig);

    for i = 1:numel(thresholdValues)
        t = thresholdValues(i);
        fig = figure('Visible', 'off', 'Name', sprintf('%s rho > %.1f', caseName, t));
        hold on; axis on; daspect([1 1 1]); view(45, 35);
        selected = rho > t;
        mdl.fe.plotSolidSelected(mdl.mesh.nodes, selected, [0.15 0.15 0.15]);
        mdl.fe.plotSolidSelected(mdl.mesh.nodes, const_elems, [0.70 0.70 0.70]);
        xlabel('x'); ylabel('y'); zlabel('z');
        title(sprintf('%s, \\rho > %.1f', strrep(caseName, '_', '\_'), t));
        saveas(fig, fullfile(caseDir, sprintf('threshold_rho_gt_%02d.png', round(10 * t))));
        close(fig);
    end
end

function plotElementDensityField(mdl, rho, const_elems)
    hold on; axis on; daspect([1 1 1]); view(45, 35);
    colormap(parula); colorbar; caxis([0 1]);

    elems = mdl.mesh.elems;
    facePattern = mdl.fe.sf.fcontours';
    faces = zeros(size(elems, 1) * size(facePattern, 1), size(facePattern, 2));
    faceColor = zeros(size(faces, 1), 1);

    row = 1;
    for e = 1:size(elems, 1)
        ef = reshape(elems(e, facePattern(:)), size(facePattern));
        n = size(ef, 1);
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

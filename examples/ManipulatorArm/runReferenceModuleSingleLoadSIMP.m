% runReferenceModuleSingleLoadSIMP
% Stage 2B: single-load SIMP compliance topology optimization on
% ReferenceModuleSolidModel. This intentionally does not implement
% multi-load p-norm aggregation or linked a/b sensitivity.

clear; close all; clc;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));
clear classes;
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');

rng(2, 'twister');

% ----- Model parameters --------------------------------------------------
E         = 2.0e9;
nu        = 0.35;
R         = 0.14;
r         = 0.08;
alpha     = 22.5;
segmentLength = 0.25;
res       = 15; %#ok<NASGU>  % kept for consistency with coupledExamples.m
res_thickness = 4;

h         = segmentLength;
alpha_deg = alpha;
res_th    = res_thickness;

% ----- Optimization parameters ------------------------------------------
loadCase      = struct('My', 1.0);
loadTol       = 0.02;
volFracTarget = 0.50;
penal         = 3.0;
Rfilter       = 1.25 * (R - r);
maxIterations = 120;
minIterations = 20;
changeTol     = 1.0e-3;
objectiveTol  = 5.0e-3;
volumeTol     = 5.0e-3;
objectiveScale = [];
moveLimit     = 0.05;
minMoveLimit  = 0.003;
moveDecay     = 0.95;
mmaDamping    = 0.25;
fdElemCount   = 8;
fdStep        = 1.0e-5;
gradTol       = 1.0e-3;

% ----- Build and validate the single pure-bending load -------------------
fprintf('Building ReferenceModuleSolidModel...\n');
mdl = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th);
mdl.analysis.isConst = true;

fprintf('  Nodes : %d\n', size(mdl.mesh.nodes, 1));
fprintf('  Elems : %d\n', size(mdl.mesh.elems, 1));
fprintf('  Load  : My = %.6g\n', loadCase.My);

mdl.applyLoadCase(loadCase);
loadValidation = mdl.validateLoadApplication(loadCase, loadTol);
assert(loadValidation.passed, 'Pure bending load validation failed.');

fixedNodeIds = find(mdl.fixedFaceSelector.select(mdl.mesh.nodes));
const_elems = find(any(ismember(mdl.mesh.elems, [mdl.loaded_node_ids(:); fixedNodeIds(:)]), 2));
const_elems = unique(const_elems(:));
fprintf('  Const elems on loaded/fixed faces: %d\n', numel(const_elems));

% ----- Create optimizer and lock the loaded/fixed face elements ----------
resCirc = round(2*pi*R / (R-r) * res_th);
Rfilter = R * 3 * pi / resCirc;
fprintf('  Rfilter from coupled model resolution: %.6g\n', Rfilter);

topOpt = SIMP_MMA_TopologyOptimizationElasticCompliance( ...
    Rfilter, mdl.analysis, penal, volFracTarget, true);
topOpt.setConstElems(const_elems);
free_elems = setdiff((1:topOpt.totalFENumber)', const_elems);
freeInitialRho = (volFracTarget * topOpt.totalFENumber - numel(const_elems)) / numel(free_elems);
assert(freeInitialRho > topOpt.xmin(1), ...
    'Requested volume fraction is infeasible with %d constant elements.', numel(const_elems));
fprintf('  Initial free-element rho: %.6f\n', freeInitialRho);
topOpt.x(:) = freeInitialRho;
topOpt.x(const_elems) = 1.0;
topOpt.xmin(const_elems) = 1.0 - 1.0e-9;
topOpt.xmax(const_elems) = 1.0;
topOpt.xold1 = topOpt.x;
topOpt.xold2 = topOpt.x;
topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
objectiveScale = 1.0 / max(abs(topOpt.FobjValue), eps);
fprintf('  MMA objective scale: %.6e\n', objectiveScale);

% ----- Finite-difference sensitivity test --------------------------------
fd = finiteDifferenceSensitivityTest(topOpt, const_elems, fdElemCount, fdStep);
fprintf('\nFinite-difference sensitivity check\n');
fprintf('  Tested elements      : %s\n', mat2str(fd.elemIds(:)'));
fprintf('  Max relative error   : %.3e\n', fd.maxRelativeError);
assert(fd.maxRelativeError < gradTol, ...
    'Finite-difference gradient error %.3e exceeds %.3e.', fd.maxRelativeError, gradTol);

% ----- Run MMA with explicit convergence history -------------------------
fprintf('\nRunning SIMP MMA compliance optimization...\n');
history = initializeHistory(maxIterations);
topOpt.iteration = 1;
topOpt.resetAnalysis();

while topOpt.iteration <= maxIterations
    previousRho = topOpt.x;
    updateDesignScaledMMA(topOpt, objectiveScale);
    currentMoveLimit = max(minMoveLimit, moveLimit * moveDecay^(topOpt.iteration - 1));
    candidateRho = topOpt.x;
    candidateRho(free_elems) = min(max(candidateRho(free_elems), ...
        previousRho(free_elems) - currentMoveLimit), previousRho(free_elems) + currentMoveLimit);
    candidateRho(const_elems) = 1.0;
    candidateRho = previousRho + mmaDamping * (candidateRho - previousRho);
    candidateRho(const_elems) = 1.0;
    topOpt.x = enforceVolumeFraction(candidateRho, free_elems, const_elems, ...
        volFracTarget, topOpt.xmin, topOpt.xmax);
    topOpt.change = max(abs(topOpt.x - previousRho));
    topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
    topOpt.computeConstraintsAndGradient(topOpt.x);
    topOpt.plotFrame();
    topOpt.printIterationInfo();

    k = topOpt.iteration;
    history.iteration(k, 1) = k;
    history.compliance(k, 1) = topOpt.FobjValue;
    history.volumeFraction(k, 1) = mean(topOpt.x);
    history.constraint(k, 1) = topOpt.constrValues;
    history.change(k, 1) = topOpt.change;

    topOpt.iteration = topOpt.iteration + 1;

    if topOpt.iteration > minIterations && ...
            (topOpt.change < changeTol || objectiveHistoryConverged(history, objectiveTol)) && ...
            abs(topOpt.constrValues) < volumeTol
        break;
    end
end

history = trimHistory(history);
topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
topOpt.computeConstraintsAndGradient(topOpt.x);

finalRho = topOpt.x;
finalCompliance = topOpt.FobjValue;
finalVolumeFraction = mean(finalRho);
finalConstraint = topOpt.constrValues;
mmaConverged = topOpt.change < changeTol || objectiveHistoryConverged(history, objectiveTol);
volumeActive = abs(finalConstraint) < volumeTol;
freeRho = finalRho(free_elems);
freeRhoStd = std(freeRho);
freeRhoRange = max(freeRho) - min(freeRho);
nonuniformTopology = freeRhoStd > 0.03 && freeRhoRange > 0.15;

fprintf('\nStage 2B checks\n');
fprintf('  Load validation passed       : %d\n', loadValidation.passed);
fprintf('  FD max relative error        : %.3e\n', fd.maxRelativeError);
fprintf('  MMA converged                : %d (change %.3e)\n', mmaConverged, topOpt.change);
fprintf('  Volume constraint active     : %d (g %.3e, vf %.4f)\n', ...
    volumeActive, finalConstraint, finalVolumeFraction);
fprintf('  Nonuniform topology          : %d (free rho std %.3e, range %.3e)\n', ...
    nonuniformTopology, freeRhoStd, freeRhoRange);

% ----- Save results ------------------------------------------------------
resultDir = fullfile(scriptDir, 'results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end
resultFile = fullfile(resultDir, 'referenceModuleSingleLoadSIMP.mat');
save(resultFile, ...
    'history', 'finalRho', 'finalVolumeFraction', 'finalCompliance', ...
    'finalConstraint', 'fd', 'loadValidation', 'const_elems', ...
    'free_elems', 'freeRhoStd', 'freeRhoRange', ...
    'volFracTarget', 'penal', 'Rfilter', 'resCirc', 'objectiveScale', ...
    'moveLimit', 'minMoveLimit', 'moveDecay', 'mmaDamping', 'loadCase', ...
    'E', 'nu', 'r', 'R', 'h', 'alpha_deg', 'res_th');
fprintf('Saved results to %s\n', resultFile);

% ----- Plot optimized density field --------------------------------------
figure('Name', 'Reference module SIMP density');
plotElementDensityField(mdl, finalRho, const_elems);
title(sprintf('Single-load SIMP density, My=1, vf=%.3f, C=%.4e', ...
    finalVolumeFraction, finalCompliance));

figure('Name', 'Reference module SIMP thresholded topology');
hold on; axis on; daspect([1 1 1]); view(45, 35);
sortedFreeRho = sort(freeRho);
denseFreeThreshold = sortedFreeRho(max(1, ceil(0.70 * numel(sortedFreeRho))));
denseFreeElems = false(size(finalRho));
denseFreeElems(free_elems) = finalRho(free_elems) >= denseFreeThreshold;
mdl.fe.plotSolidSelected(mdl.mesh.nodes, denseFreeElems, [0.15 0.15 0.15]);
mdl.fe.plotSolidSelected(mdl.mesh.nodes, const_elems, [0.70 0.70 0.70]);
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('Densest free elements, rho >= %.3f', denseFreeThreshold));

figure('Name', 'SIMP convergence history');
tiledlayout(2, 1);
nexttile;
plot(history.iteration, history.compliance, '-o', 'LineWidth', 1.0);
grid on; xlabel('Iteration'); ylabel('Compliance');
nexttile;
plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.0);
hold on; yline(volFracTarget, '--');
grid on; xlabel('Iteration'); ylabel('Volume fraction');

assert(mmaConverged, 'MMA did not converge within %d iterations.', maxIterations);
assert(volumeActive, 'Volume constraint is not active enough: g = %.3e.', finalConstraint);
assert(nonuniformTopology, 'Optimized density field is too uniform.');

% =========================================================================
function history = initializeHistory(n)
    history.iteration = nan(n, 1);
    history.compliance = nan(n, 1);
    history.volumeFraction = nan(n, 1);
    history.constraint = nan(n, 1);
    history.change = nan(n, 1);
end

function history = trimHistory(history)
    keep = ~isnan(history.iteration);
    fields = fieldnames(history);
    for i = 1:numel(fields)
        history.(fields{i}) = history.(fields{i})(keep);
    end
end

function tf = objectiveHistoryConverged(history, tol)
    compliance = history.compliance(~isnan(history.compliance));
    nWindow = min(10, numel(compliance));
    if nWindow < 3
        tf = false;
        return;
    end
    c = compliance(end-nWindow+1:end);
    tf = (max(c) - min(c)) / max(abs(c(end)), eps) < tol;
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

function fd = finiteDifferenceSensitivityTest(topOpt, const_elems, nTests, h)
    x0 = topOpt.x;
    freeIds = setdiff((1:numel(x0))', const_elems(:));
    sampleCount = min(nTests, numel(freeIds));
    elemIds = freeIds(randperm(numel(freeIds), sampleCount));

    topOpt.x = x0;
    topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
    analyticGrad = topOpt.gradFobjValue;

    fdGrad = zeros(sampleCount, 1);
    relErr = zeros(sampleCount, 1);

    for i = 1:sampleCount
        eid = elemIds(i);
        step = min([h, 0.49 * (x0(eid) - topOpt.xmin(eid)), 0.49 * (topOpt.xmax(eid) - x0(eid))]);
        if step <= 0
            error('Element %d has no room for finite-difference perturbation.', eid);
        end

        xp = x0;
        xm = x0;
        xp(eid) = xp(eid) + step;
        xm(eid) = xm(eid) - step;

        topOpt.x = xp;
        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        cp = topOpt.FobjValue;

        topOpt.x = xm;
        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        cm = topOpt.FobjValue;

        fdGrad(i) = (cp - cm) / (2 * step);
        relErr(i) = abs(fdGrad(i) - analyticGrad(eid)) / ...
            max([abs(fdGrad(i)), abs(analyticGrad(eid)), eps]);
    end

    topOpt.x = x0;
    topOpt.computeObjectiveFunctonWithGradient(topOpt.x);

    fd.elemIds = elemIds;
    fd.analyticGradient = analyticGrad(elemIds);
    fd.finiteDifferenceGradient = fdGrad;
    fd.relativeError = relErr;
    fd.maxRelativeError = max(relErr);
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

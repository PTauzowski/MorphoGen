% runReferenceModuleMultiLoadSIMP
% Stage 2C: multi-load p-norm SIMP compliance optimization on
% ReferenceModuleSolidModel using synthetic pure load cases only.

clear; close all; clc;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));
clear classes;
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..'); %#ok<NASGU>

rng(3, 'twister');

% ----- Model parameters from coupledExamples.m ---------------------------
E = 2.0e9;
nu = 0.35;
R = 0.14;
r = 0.08;
alpha = 22.5;
segmentLength = 0.25;
res = 15; %#ok<NASGU>
res_thickness = 4;

h = segmentLength;
alpha_deg = alpha;
res_th = res_thickness;

% ----- Optimization parameters ------------------------------------------
loadTol = 0.02;
volFracTarget = 0.35;
penal = 4.0;
pAgg = 6.0;
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
fdElemCount = 8;
fdStep = 1.0e-5;
gradTol = 1.0e-3;
dominanceTol = 0.90;

loadNames = ["My=1"; "Mz=1"; "Ms=1"];
loadCases = {
    struct('My', 1.0);
    struct('Mz', 1.0);
    struct('Ms', 1.0);
    % Add shear cases after the bending/torsion benchmark is understood.
    % struct('Ty', 1.0);
    % struct('Tz', 1.0);
};
weights = ones(numel(loadCases), 1) / numel(loadCases);

% ----- Build and validate all synthetic loads ----------------------------
fprintf('Building ReferenceModuleSolidModel...\n');
mdl = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th);
mdl.analysis.isConst = true;

fprintf('  Nodes : %d\n', size(mdl.mesh.nodes, 1));
fprintf('  Elems : %d\n', size(mdl.mesh.elems, 1));

loadValidation = cell(numel(loadCases), 1);
allLoadsPassed = true;
for k = 1:numel(loadCases)
    fprintf('\n====== Load case: %s ======\n', char(loadNames(k)));
    mdl.applyLoadCase(loadCases{k});
    loadValidation{k} = mdl.validateLoadApplication(loadCases{k}, loadTol);
    allLoadsPassed = allLoadsPassed && loadValidation{k}.passed;
end
assert(allLoadsPassed, 'At least one synthetic load case failed validation.');

fixedNodeIds = find(mdl.fixedFaceSelector.select(mdl.mesh.nodes));
const_elems = find(any(ismember(mdl.mesh.elems, [mdl.loaded_node_ids(:); fixedNodeIds(:)]), 2));
const_elems = unique(const_elems(:));
fprintf('\nConst elems on loaded/fixed faces: %d\n', numel(const_elems));

resCirc = round(2*pi*R / (R-r) * res_th);
Rfilter = R * 3 * pi / resCirc;
Rfilter = 0.60 * Rfilter;
fprintf('Rfilter from coupled model resolution: %.6g\n', Rfilter);

% ----- Create optimizer and initialize normalized multi-load objective ----
topOpt = SIMP_MMA_ReferenceModuleMultiLoadCompliance( ...
    Rfilter, mdl, loadCases, weights, pAgg, penal, volFracTarget, true, useComplianceNormalization);
topOpt.setConstElems(const_elems);

free_elems = setdiff((1:topOpt.totalFENumber)', const_elems);
freeInitialRho = (volFracTarget * topOpt.totalFENumber - numel(const_elems)) / numel(free_elems);
assert(freeInitialRho > topOpt.xmin(1), ...
    'Requested volume fraction is infeasible with %d constant elements.', numel(const_elems));
fprintf('Initial free-element rho: %.6f\n', freeInitialRho);

topOpt.x(:) = freeInitialRho;
topOpt.x(const_elems) = 1.0;
topOpt.xmin(const_elems) = 1.0 - 1.0e-9;
topOpt.xmax(const_elems) = 1.0;
topOpt.xold1 = topOpt.x;
topOpt.xold2 = topOpt.x;
topOpt.initializeComplianceNormalization(topOpt.x);
topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
objectiveScale = 1.0 / max(abs(topOpt.FobjValue), eps);

fprintf('Initial compliance normalization:\n');
for k = 1:numel(loadCases)
    fprintf('  %-6s C0 = %.6e\n', char(loadNames(k)), topOpt.complianceNormalization(k));
end
fprintf('MMA objective scale: %.6e\n', objectiveScale);

% ----- Finite-difference sensitivity test for aggregate objective --------
fd = finiteDifferenceSensitivityTest(topOpt, const_elems, fdElemCount, fdStep);
fprintf('\nAggregated finite-difference sensitivity check\n');
fprintf('  Tested elements      : %s\n', mat2str(fd.elemIds(:)'));
fprintf('  Max relative error   : %.3e\n', fd.maxRelativeError);
assert(fd.maxRelativeError < gradTol, ...
    'Aggregated finite-difference gradient error %.3e exceeds %.3e.', fd.maxRelativeError, gradTol);

% ----- Run MMA with explicit convergence history -------------------------
fprintf('\nRunning multi-load p-norm SIMP MMA optimization...\n');
history = initializeHistory(maxIterations, numel(loadCases));
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
    history.objective(k, 1) = topOpt.FobjValue;
    history.compliance(k, :) = topOpt.complianceValues(:)';
    history.normalizedCompliance(k, :) = topOpt.normalizedComplianceValues(:)';
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
finalObjective = topOpt.FobjValue;
finalCompliance = topOpt.complianceValues;
finalNormalizedCompliance = topOpt.normalizedComplianceValues;
finalVolumeFraction = mean(finalRho);
finalConstraint = topOpt.constrValues;
mmaConverged = topOpt.change < changeTol || objectiveHistoryConverged(history, objectiveTol);
volumeActive = abs(finalConstraint) < volumeTol;
freeRho = finalRho(free_elems);
freeRhoStd = std(freeRho);
freeRhoRange = max(freeRho) - min(freeRho);
nonuniformTopology = freeRhoStd > 0.03 && freeRhoRange > 0.15;
loadContributions = weights .* (finalNormalizedCompliance .^ pAgg);
loadContributions = loadContributions / sum(loadContributions);
noDominatingLoad = max(loadContributions) < dominanceTol;

fprintf('\nStage 2C checks\n');
fprintf('  Every load validation passed : %d\n', allLoadsPassed);
fprintf('  Aggregate FD max rel error   : %.3e\n', fd.maxRelativeError);
fprintf('  MMA converged                : %d (change %.3e)\n', mmaConverged, topOpt.change);
fprintf('  Volume constraint active     : %d (g %.3e, vf %.4f)\n', ...
    volumeActive, finalConstraint, finalVolumeFraction);
fprintf('  Nonuniform topology          : %d (free rho std %.3e, range %.3e)\n', ...
    nonuniformTopology, freeRhoStd, freeRhoRange);
fprintf('  No normalized load dominance : %d (max contribution %.3f)\n', ...
    noDominatingLoad, max(loadContributions));
for k = 1:numel(loadCases)
    fprintf('    %-6s C=%.6e  C/C0=%.6f  contribution=%.3f\n', ...
        char(loadNames(k)), finalCompliance(k), finalNormalizedCompliance(k), loadContributions(k));
end

% ----- Save results ------------------------------------------------------
resultDir = fullfile(scriptDir, 'results');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end
resultFile = fullfile(resultDir, 'referenceModuleMultiLoadSIMP.mat');
save(resultFile, ...
    'history', 'finalRho', 'finalObjective', 'finalCompliance', ...
    'finalNormalizedCompliance', 'finalVolumeFraction', 'finalConstraint', ...
    'fd', 'loadValidation', 'loadCases', 'loadNames', 'weights', ...
    'loadContributions', 'const_elems', 'free_elems', 'freeRhoStd', ...
    'freeRhoRange', 'volFracTarget', 'penal', 'pAgg', 'Rfilter', ...
    'resCirc', 'objectiveScale', 'moveLimit', 'minMoveLimit', ...
    'moveDecay', 'mmaDamping', 'useComplianceNormalization', ...
    'E', 'nu', 'r', 'R', 'h', 'alpha_deg', 'res_th');
fprintf('Saved results to %s\n', resultFile);

% ----- Plot optimized density field and compliance history ---------------
figure('Name', 'Reference module multi-load SIMP density');
plotElementDensityField(mdl, finalRho, const_elems);
title(sprintf('Multi-load p-norm SIMP density, vf=%.3f, J=%.4e', ...
    finalVolumeFraction, finalObjective));

figure('Name', 'Reference module multi-load thresholded topology');
hold on; axis on; daspect([1 1 1]); view(45, 35);
sortedFreeRho = sort(freeRho);
denseFreeThreshold = sortedFreeRho(max(1, ceil(0.80 * numel(sortedFreeRho))));
denseFreeElems = false(size(finalRho));
denseFreeElems(free_elems) = finalRho(free_elems) >= denseFreeThreshold;
mdl.fe.plotSolidSelected(mdl.mesh.nodes, denseFreeElems, [0.15 0.15 0.15]);
mdl.fe.plotSolidSelected(mdl.mesh.nodes, const_elems, [0.70 0.70 0.70]);
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('Densest free elements, rho >= %.3f', denseFreeThreshold));

figure('Name', 'Multi-load SIMP convergence history');
tiledlayout(3, 1);
nexttile;
plot(history.iteration, history.objective, '-o', 'LineWidth', 1.0);
grid on; xlabel('Iteration'); ylabel('p-norm J');
nexttile;
plot(history.iteration, history.normalizedCompliance, 'LineWidth', 1.0);
grid on; xlabel('Iteration'); ylabel('C_k / C_{k0}');
legend(loadNames, 'Location', 'best');
nexttile;
plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.0);
hold on; yline(volFracTarget, '--');
grid on; xlabel('Iteration'); ylabel('Volume fraction');

assert(allLoadsPassed, 'At least one load validation failed.');
assert(fd.maxRelativeError < gradTol, ...
    'Aggregated finite-difference gradient error %.3e exceeds %.3e.', fd.maxRelativeError, gradTol);
assert(mmaConverged, 'MMA did not converge within %d iterations.', maxIterations);
assert(volumeActive, 'Volume constraint is not active enough: g = %.3e.', finalConstraint);
assert(nonuniformTopology, 'Optimized density field is too uniform.');
assert(noDominatingLoad, 'A normalized load contribution dominates the p-norm aggregate.');

% =========================================================================
function history = initializeHistory(nIter, nLoads)
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
        history.(fields{i}) = history.(fields{i})(keep, :);
    end
end

function tf = objectiveHistoryConverged(history, tol)
    objective = history.objective(~isnan(history.objective));
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
        jp = topOpt.FobjValue;

        topOpt.x = xm;
        topOpt.computeObjectiveFunctonWithGradient(topOpt.x);
        jm = topOpt.FobjValue;

        fdGrad(i) = (jp - jm) / (2 * step);
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

% testA_fullArmUnlinkedSIMP
% Test A: full-arm unlinked multi-configuration SIMP diagnostic.
%
% Diagnostic purpose:
%   Optimize the full 3D Arm-Z arm directly with one independent density
%   variable per full-arm solid element. This intentionally does not enforce
%   repeated-module manufacturability or segmentToArm linking. It is an
%   upper-bound diagnostic for what an unlinked full-arm design can do.
%
% Configurations:
%   1. Max M_z bending
%   2. Max M_s torsion
%   3. Max T_y shear
%
% Objective:
%   J = (sum_k w_k * (C_k/C0_k)^pAgg)^(1/pAgg)
%
% where C_k = f_k' * u_k and C0_k is the compliance of the initial
% VolFrac density field.

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(11, 'twister');

%% ---- Geometry and optimization parameters ------------------------------
E = 2.0e9;
nu = 0.35;
R = 0.14;
r = 0.08;
h_seg = 0.25;
alpha = 22.5;
res = 15;
res_th = 4;
Pz = 100;
ShapeFn = ShapeFunctionL8();

VolFrac = 0.40;
penal = 3.0;
pAgg = 4.0;
maxIter = 30;
xminValue = 0.01;

fdElemCount = 8;
fdStep = 1.0e-5;
fdRelTol = 5.0e-3;   % p-norm aggregation + 48k elems: FD noise ~ 4e-3 is expected

moveLimit = 0.05;
minMoveLimit = 0.003;
moveDecay = 0.95;
mmaDamping = 0.25;
objectiveTol = 5.0e-3;

sampleMinN_smooth=[0 180 180 180 180 180 180]; %  0  164.4607  177.0448  215.1802  214.3430  186.1769  240.6389
sampleMaxMz_smooth=[0 0 0 180 180 180 180];     %  0  343.7190   54.9144  125.4541  163.9858  169.3933  110.6999
sampleMaxTy_smooth=[ 0  0  180 0 180  180  180]; %  0    0.0253  246.2603  287.9922  340.8601  106.8641  112.8238
sampleMaxMs_smooth=[ 0  45 45 45  270  180 180]; %  0 0   33.5628   56.5082   36.2258  278.0525  195.3981  122.2441

configs = {
    struct('name', 'min_bending', 'label', 'Min M_z', 'betas', -[0 0 0 180 180 180 180]);
    struct('name', 'min_torsion', 'label', 'Min M_s', 'betas', -[0 45 45 45 270 180 180]);
    struct('name', 'min_shear',   'label', 'Min T_y', 'betas', -[0 0 180 0 180 180 180]);
    struct('name', 'max_bending', 'label', 'Max M_z', 'betas', [0 0 0 180 180 180 180]);
    struct('name', 'max_torsion', 'label', 'Max M_s', 'betas', [0 45 45 45 270 180 180]);
    struct('name', 'max_shear',   'label', 'Max T_y', 'betas', [0 0 180 0 180 180 180]);
};
nConfigs = numel(configs);
weights = ones(nConfigs, 1) / nConfigs;

resultRoot = fullfile(scriptDir, 'results', 'testA_fullArmUnlinkedSIMP');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test A full-arm unlinked SIMP diagnostic\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf('  VolFrac=%.3f, penal=%.2f, pAgg=%.2f, maxIter=%d, xmin=%.3f\n', ...
    VolFrac, penal, pAgg, maxIter, xminValue);

%% ---- Build full-arm configurations -------------------------------------
models = cell(nConfigs, 1);
analyses = cell(nConfigs, 1);
setupRows = cell(nConfigs, 1);

referenceElems = [];
referenceElemCount = [];
referenceDofs = [];
referenceTaskDim = [];

for k = 1:nConfigs
    cfg = configs{k};
    fprintf('\nBuilding configuration %d/%d: %s, betas=%s\n', ...
        k, nConfigs, cfg.label, mat2str(cfg.betas));

    model = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
        cfg.betas, ShapeFn, true, Pz);
    analysis = model.analysis;
    nElems = analysis.getTotalElemsNumber();
    taskDim = analysis.getTaskDim();
    nSupports = nnz(analysis.supports);
    nLoadedDofs = nnz(analysis.Pnodal);

    if k == 1
        referenceElems = model.mesh.elems;
        referenceElemCount = nElems;
        referenceDofs = analysis.ndofs;
        referenceTaskDim = taskDim;
        fprintf('  Reference element count: %d\n', referenceElemCount);
        fprintf('  Reference task DOFs    : %d\n', referenceTaskDim);
    else
        assert(nElems == referenceElemCount, ...
            'Element-count mismatch in %s: got %d, expected %d.', ...
            cfg.name, nElems, referenceElemCount);
        assert(taskDim == referenceTaskDim, ...
            'DOF-count mismatch in %s: got %d, expected %d.', ...
            cfg.name, taskDim, referenceTaskDim);
        assert(isequal(model.mesh.elems, referenceElems), ...
            'Mesh connectivity differs in configuration %s.', cfg.name);
        assert(isequal(analysis.ndofs, referenceDofs), ...
            'DOF labels/order differ in configuration %s.', cfg.name);
    end

    models{k} = model;
    analyses{k} = analysis;

    row.configName = string(cfg.name);
    row.configLabel = string(cfg.label);
    row.nNodes = size(model.mesh.nodes, 1);
    row.nElems = nElems;
    row.nTaskDofs = taskDim;
    row.nSupportedDofs = nSupports;
    row.nLoadedDofsBeforeSolve = nLoadedDofs;
    row.sameConnectivityAsFirst = true;
    row.sameDofsAsFirst = true;
    setupRows{k, 1} = row;

    fprintf('  nodes=%d, elems=%d, taskDOFs=%d, supportedDOFs=%d, loadedDOFs=%d\n', ...
        row.nNodes, row.nElems, row.nTaskDofs, row.nSupportedDofs, row.nLoadedDofsBeforeSolve);
end

setupTable = struct2table(vertcat(setupRows{:}));
writetable(setupTable, fullfile(resultRoot, 'configuration_setup.csv'));

nDesign = referenceElemCount;
xmin = xminValue * ones(nDesign, 1);
xmax = ones(nDesign, 1);
x = VolFrac * ones(nDesign, 1);
x = enforceVolumeFraction(x, VolFrac, xmin, xmax);

%% ---- Compliance normalization and finite-difference gradient check -------
fprintf('\nComputing initial compliance normalization at x=VolFrac...\n');
[J0, grad0, C0, Cinit] = evaluateObjectiveAndGradient(analyses, x, penal, pAgg, weights, []);
fprintf('  Initial J = %.8e\n', J0);
for k = 1:nConfigs
    fprintf('    %-12s C0 = %.8e\n', configs{k}.name, C0(k));
end

fprintf('\nFinite-difference gradient test (%d random elements, h=%.1e)...\n', ...
    fdElemCount, fdStep);
fd = finiteDifferenceGradientTest(analyses, x, penal, pAgg, weights, C0, ...
    grad0, fdElemCount, fdStep, xmin, xmax);
fprintf('  FD max relative error = %.3e\n', fd.maxRelativeError);
assert(fd.maxRelativeError < fdRelTol, ...
    'FD gradient check failed: max relative error %.3e exceeds %.3e.', ...
    fd.maxRelativeError, fdRelTol);

fdTable = struct2table(fd.rows);
writetable(fdTable, fullfile(resultRoot, 'finite_difference_gradient.csv'));

%% ---- MMA optimization ---------------------------------------------------
fprintf('\nRunning projected MMA for %d iterations...\n', maxIter);

xHistory = zeros(nDesign, maxIter + 1);
xHistory(:, 1) = x;
JHistory = nan(maxIter + 1, 1);
JHistory(1) = J0;
CHistory = nan(maxIter + 1, nConfigs);
CHistory(1, :) = Cinit(:)';
volHistory = nan(maxIter + 1, 1);
volHistory(1) = mean(x);
changeHistory = nan(maxIter + 1, 1);
changeHistory(1) = 0;

m = 1;
n = nDesign;
xold1 = x;
xold2 = x;
low = zeros(n, 1);
upp = ones(n, 1);
a0 = 1;
a = 0;
c = 1000;
d = 0;
objectiveScale = 1.0 / max(abs(J0), eps);

for iter = 1:maxIter
    [J, gradJ, ~, ~] = evaluateObjectiveAndGradient(analyses, x, penal, pAgg, weights, C0);
    constr = sum(x) / (VolFrac * nDesign) - 1.0;
    gradConstr = ones(1, nDesign) / (VolFrac * nDesign);

    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = mmasub2( ...
        m, n, iter, x, xmin, xmax, xold1, xold2, ...
        objectiveScale * J, objectiveScale * gradJ, 0 * gradJ, ...
        constr, gradConstr, 0 * gradConstr, ...
        low, upp, a0, a, c, d);

    if iter > 1
        xold2 = xold1;
    end
    xold1 = x;

    currentMoveLimit = max(minMoveLimit, moveLimit * moveDecay^(iter - 1));
    xCandidate = min(max(xmma, x - currentMoveLimit), x + currentMoveLimit);
    xCandidate = x + mmaDamping * (xCandidate - x);
    xCandidate = enforceVolumeFraction(xCandidate, VolFrac, xmin, xmax);

    change = max(abs(xCandidate - x));
    x = xCandidate;

    [Jnew, Cnew] = evaluateObjectiveOnly(analyses, x, penal, pAgg, weights, C0);
    xHistory(:, iter + 1) = x;
    JHistory(iter + 1) = Jnew;
    CHistory(iter + 1, :) = Cnew(:)';
    volHistory(iter + 1) = mean(x);
    changeHistory(iter + 1) = change;

    fprintf('%4d  J=%12.6e  dJ/J0=% .3e  vf=%.4f  change=%.3e', ...
        iter, Jnew, (Jnew - JHistory(iter)) / max(abs(JHistory(iter)), eps), ...
        volHistory(iter + 1), change);
    for k = 1:nConfigs
        fprintf('  C_%s=%.3e', configs{k}.name, Cnew(k));
    end
    fprintf('\n');
end

finalX = x;
finalJ = JHistory(maxIter + 1);
finalC = CHistory(maxIter + 1, :)';
finalVolumeFraction = mean(finalX);
finalChange = changeHistory(maxIter + 1);
objectiveWindowConverged = objectiveHistoryConverged(JHistory, objectiveTol);
objectiveDecreased = finalJ < JHistory(1);
mmaAcceptable = objectiveWindowConverged || objectiveDecreased;
volumeActive = abs(finalVolumeFraction - VolFrac) < 5.0e-3;
topologyNonuniform = std(finalX) > 0.03 && (max(finalX) - min(finalX)) > 0.15;

%% ---- Location-variation diagnostics ------------------------------------
locationStats = computeLocationDensityStats(models{1}, finalX);
locationDensityRange = max(locationStats.meanDensity) - min(locationStats.meanDensity);
locationDensityStd = std(locationStats.meanDensity);
topologyDiffersByArmLocation = locationDensityRange > 0.03;

writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

fprintf('\nFinal checks\n');
fprintf('  Objective decreased          : %d (J0=%.6e, Jf=%.6e)\n', ...
    objectiveDecreased, JHistory(1), finalJ);
fprintf('  Objective window converged   : %d\n', objectiveWindowConverged);
fprintf('  MMA acceptable               : %d\n', mmaAcceptable);
fprintf('  Volume constraint active     : %d (vf=%.6f)\n', volumeActive, finalVolumeFraction);
fprintf('  Topology nonuniform          : %d (std=%.4f, range=%.4f)\n', ...
    topologyNonuniform, std(finalX), max(finalX) - min(finalX));
fprintf('  Differs by arm location      : %d (half-seg mean range=%.4f)\n', ...
    topologyDiffersByArmLocation, locationDensityRange);

%% ---- Save histories and plots ------------------------------------------
history.iteration = (0:maxIter)';
history.J = JHistory;
history.C = CHistory;
history.volumeFraction = volHistory;
history.change = changeHistory;
history.x = xHistory;
history.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
history.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));
history.weights = weights;
history.C0 = C0;

save(fullfile(resultRoot, 'result.mat'), ...
    'finalX', 'finalJ', 'finalC', 'finalVolumeFraction', 'finalChange', ...
    'history', 'fd', 'configs', 'VolFrac', 'penal', 'pAgg', 'maxIter', ...
    'xminValue', 'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', ...
    'Pz', 'locationStats', 'locationDensityRange', 'locationDensityStd', ...
    'topologyDiffersByArmLocation', 'objectiveDecreased', ...
    'objectiveWindowConverged', 'mmaAcceptable', 'volumeActive', ...
    'topologyNonuniform');

saveHistoryCsv(resultRoot, history, configs);
saveSummaryCsv(resultRoot, fd, history, finalX, finalC, VolFrac, ...
    objectiveDecreased, objectiveWindowConverged, mmaAcceptable, volumeActive, ...
    topologyNonuniform, topologyDiffersByArmLocation, ...
    locationDensityRange, locationDensityStd);

plotFinalTopology(models{1}, finalX, resultRoot);
plotHistory(history, configs, resultRoot);
plotLocationDensity(locationStats, resultRoot);

fprintf('\nSaved Test A outputs to %s\n', resultRoot);

assert(mmaAcceptable, ...
    'MMA did not converge by objective window and did not decrease J over %d iterations.', maxIter);
assert(volumeActive, ...
    'Volume constraint is not active enough: vf=%.6f, target=%.6f.', finalVolumeFraction, VolFrac);
assert(topologyNonuniform, ...
    'Final topology is too uniform: std=%.4f, range=%.4f.', ...
    std(finalX), max(finalX) - min(finalX));

%% =========================================================================
function [J, gradJ, C0, C] = evaluateObjectiveAndGradient(analyses, x, penal, pAgg, weights, C0)
    nConfigs = numel(analyses);
    C = zeros(nConfigs, 1);
    dC = zeros(numel(x), nConfigs);

    for k = 1:nConfigs
        [C(k), dC(:, k)] = computeComplianceAndGradient(analyses{k}, x, penal);
    end

    if isempty(C0)
        C0 = max(C, eps);
    end

    Cagg = C ./ C0;
    weightedSum = sum(weights(:) .* (Cagg .^ pAgg));
    J = weightedSum ^ (1.0 / pAgg);

    gradJ = zeros(numel(x), 1);
    if J > eps
        for k = 1:nConfigs
            gradJ = gradJ + weights(k) * Cagg(k)^(pAgg - 1) * dC(:, k) / C0(k);
        end
        gradJ = J^(1 - pAgg) * gradJ;
    end
end

function [J, C] = evaluateObjectiveOnly(analyses, x, penal, pAgg, weights, C0)
    nConfigs = numel(analyses);
    C = zeros(nConfigs, 1);

    for k = 1:nConfigs
        C(k) = computeComplianceOnly(analyses{k}, x, penal);
    end

    assert(~isempty(C0), 'Compliance normalization C0 must be provided.');
    Cagg = C ./ C0;
    weightedSum = sum(weights(:) .* (Cagg .^ pAgg));
    J = weightedSum ^ (1.0 / pAgg);
end

function [C, dC] = computeComplianceAndGradient(analysis, x, penal)
    nElemsTotal = analysis.getTotalElemsNumber();
    assert(numel(x) == nElemsTotal, ...
        'Density vector length %d does not match analysis element count %d.', ...
        numel(x), nElemsTotal);

    x = x(:);
    qfem = analysis.solveWeighted(x .^ penal);
    q = analysis.fromFEMVector(qfem(:, 1));
    P = analysis.Pfem(:, 1);
    C = P' * qfem(:, 1);

    if analysis.isConst
        stiffnessFunction = 'computeStifnessMatrixConst';
    else
        stiffnessFunction = 'computeStifnessMatrix';
    end

    dC = zeros(nElemsTotal, 1);
    elemOffset = 0;
    xOnes = ones(nElemsTotal, 1);
    elemIndices = analysis.getElemIndices();

    for i = 1:numel(analysis.felems)
        fe = analysis.felems{i};
        elemIds = elemIndices{i};
        nelems = size(fe.elems, 1);
        nnodes = size(fe.elems, 2);
        ndofs = size(fe.ndofs, 2);
        dim = nnodes * ndofs;
        K0 = reshape(fe.(stiffnessFunction)(analysis.mesh.nodes, xOnes(elemIds)), ...
            dim, dim, nelems);
        qelems = fe.createElemSolutionVectors(q);

        for e = 1:nelems
            globalElem = elemOffset + e;
            elemEnergy = qelems(:, e)' * K0(:, :, e) * qelems(:, e);
            dC(globalElem) = -penal * x(globalElem)^(penal - 1) * elemEnergy;
        end
        elemOffset = elemOffset + nelems;
    end
end

function C = computeComplianceOnly(analysis, x, penal)
    nElemsTotal = analysis.getTotalElemsNumber();
    assert(numel(x) == nElemsTotal, ...
        'Density vector length %d does not match analysis element count %d.', ...
        numel(x), nElemsTotal);

    qfem = analysis.solveWeighted(x(:) .^ penal);
    P = analysis.Pfem(:, 1);
    C = P' * qfem(:, 1);
end

function fd = finiteDifferenceGradientTest(analyses, x, penal, pAgg, weights, C0, grad, nTest, h, xmin, xmax)
    n = numel(x);
    nTest = min(nTest, n);
    elemIds = randperm(n, nTest)';
    rows = repmat(struct('elemId', 0, 'analytic', 0, 'finiteDifference', 0, ...
        'relativeError', 0), nTest, 1);
    relErrors = zeros(nTest, 1);

    for i = 1:nTest
        e = elemIds(i);
        hUse = min([h, 0.49 * (xmax(e) - x(e)), 0.49 * (x(e) - xmin(e))]);
        if hUse <= 0
            hUse = h;
        end

        xp = x;
        xm = x;
        xp(e) = min(xmax(e), xp(e) + hUse);
        xm(e) = max(xmin(e), xm(e) - hUse);

        [Jp, ~] = evaluateObjectiveOnly(analyses, xp, penal, pAgg, weights, C0);
        [Jm, ~] = evaluateObjectiveOnly(analyses, xm, penal, pAgg, weights, C0);
        fdGrad = (Jp - Jm) / (xp(e) - xm(e));
        relErr = abs(fdGrad - grad(e)) / max([abs(fdGrad), abs(grad(e)), eps]);

        rows(i).elemId = e;
        rows(i).analytic = grad(e);
        rows(i).finiteDifference = fdGrad;
        rows(i).relativeError = relErr;
        relErrors(i) = relErr;

        fprintf('  elem %7d: analytic=% .6e  FD=% .6e  relErr=%.3e\n', ...
            e, grad(e), fdGrad, relErr);
    end

    fd.elemIds = elemIds;
    fd.rows = rows;
    fd.maxRelativeError = max(relErrors);
    fd.meanRelativeError = mean(relErrors);
end

function x = enforceVolumeFraction(x, VolFrac, xmin, xmax)
    x = min(max(x(:), xmin), xmax);
    target = VolFrac * numel(x);
    target = min(max(target, sum(xmin)), sum(xmax));

    lo = min(xmin - x);
    hi = max(xmax - x);
    for k = 1:80
        shift = 0.5 * (lo + hi);
        candidate = min(max(x + shift, xmin), xmax);
        if sum(candidate) < target
            lo = shift;
        else
            hi = shift;
        end
    end
    x = min(max(x + 0.5 * (lo + hi), xmin), xmax);
end

function tf = objectiveHistoryConverged(JHistory, tol)
    objective = JHistory(~isnan(JHistory));
    nWindow = min(10, numel(objective));
    if nWindow < 3
        tf = false;
        return;
    end
    j = objective(end-nWindow+1:end);
    tf = (max(j) - min(j)) / max(abs(j(end)), eps) < tol;
end

function stats = computeLocationDensityStats(model, x)
    nHalf = model.halfSegmentNelems;
    nLocations = floor(numel(x) / nHalf);
    assert(nLocations * nHalf == numel(x), ...
        'Full-arm element count %d is not divisible by half-segment count %d.', ...
        numel(x), nHalf);

    stats = repmat(struct('locationIndex', 0, 'elemStart', 0, 'elemEnd', 0, ...
        'meanDensity', 0, 'stdDensity', 0, 'rhoGt05', 0, 'rhoGt07', 0), ...
        nLocations, 1);
    for i = 1:nLocations
        ids = ((i - 1) * nHalf + 1):(i * nHalf);
        xi = x(ids);
        stats(i).locationIndex = i;
        stats(i).elemStart = ids(1);
        stats(i).elemEnd = ids(end);
        stats(i).meanDensity = mean(xi);
        stats(i).stdDensity = std(xi);
        stats(i).rhoGt05 = mean(xi > 0.5);
        stats(i).rhoGt07 = mean(xi > 0.7);
    end
end

function saveHistoryCsv(resultRoot, history, configs)
    nRows = numel(history.iteration);
    rows = repmat(struct(), nRows, 1);
    for i = 1:nRows
        rows(i).iteration = history.iteration(i);
        rows(i).J = history.J(i);
        rows(i).volumeFraction = history.volumeFraction(i);
        rows(i).change = history.change(i);
        for k = 1:numel(configs)
            fieldName = ['C_' configs{k}.name];
            rows(i).(fieldName) = history.C(i, k);
        end
    end
    writetable(struct2table(rows), fullfile(resultRoot, 'history.csv'));
end

function saveSummaryCsv(resultRoot, fd, history, finalX, finalC, VolFrac, ...
        objectiveDecreased, objectiveWindowConverged, mmaAcceptable, volumeActive, ...
        topologyNonuniform, topologyDiffersByArmLocation, ...
        locationDensityRange, locationDensityStd)
    row.status = "ok";
    row.allConfigurationsSolved = true;
    row.nElements = numel(finalX);
    row.nIterations = numel(history.iteration) - 1;
    row.fdMaxRelativeError = fd.maxRelativeError;
    row.fdPass = fd.maxRelativeError < 5.0e-3;
    row.initialJ = history.J(1);
    row.finalJ = history.J(end);
    row.relativeJChange = (history.J(end) - history.J(1)) / max(abs(history.J(1)), eps);
    row.targetVolumeFraction = VolFrac;
    row.finalVolumeFraction = mean(finalX);
    row.volumeActive = volumeActive;
    row.finalChange = history.change(end);
    row.objectiveDecreased = objectiveDecreased;
    row.objectiveWindowConverged = objectiveWindowConverged;
    row.mmaAcceptable = mmaAcceptable;
    row.topologyNonuniform = topologyNonuniform;
    row.finalDensityStd = std(finalX);
    row.finalDensityRange = max(finalX) - min(finalX);
    row.rhoGt05 = mean(finalX > 0.5);
    row.rhoGt07 = mean(finalX > 0.7);
    row.finalCMaxBending = finalC(1);
    row.finalCMaxTorsion = finalC(2);
    row.finalCMaxShear = finalC(3);
    row.topologyDiffersByArmLocation = topologyDiffersByArmLocation;
    row.locationMeanDensityRange = locationDensityRange;
    row.locationMeanDensityStd = locationDensityStd;
    writetable(struct2table(row), fullfile(resultRoot, 'summary.csv'));
end

function plotFinalTopology(model, x, resultRoot)
    fig = figure('Visible', 'off', 'Name', 'Test A final density');
    plotElementDensityField(model, x);
    title(sprintf('Test A final density, vf=%.3f', mean(x)));
    saveas(fig, fullfile(resultRoot, 'final_density.png'));
    close(fig);

    thresholds = [0.5, 0.7];
    for i = 1:numel(thresholds)
        t = thresholds(i);
        fig = figure('Visible', 'off', 'Name', sprintf('Test A rho > %.1f', t));
        hold on; axis on; daspect([1 1 1]); view(45, 35);
        selected = x > t;
        model.fe.plotSolidSelected(model.mesh.nodes, selected, [0.15 0.15 0.15]);
        xlabel('x'); ylabel('y'); zlabel('z');
        title(sprintf('Test A final topology, \\rho > %.1f', t));
        saveas(fig, fullfile(resultRoot, sprintf('final_threshold_rho_gt_%02d.png', round(10 * t))));
        close(fig);
    end
end

function plotElementDensityField(model, x)
    hold on; axis on; daspect([1 1 1]); view(45, 35);
    colormap(parula); colorbar; caxis([0 1]);

    elems = model.mesh.elems;
    facePattern = model.fe.sf.fcontours';
    faces = zeros(size(elems, 1) * size(facePattern, 1), size(facePattern, 2));
    faceColor = zeros(size(faces, 1), 1);

    row = 1;
    for e = 1:size(elems, 1)
        ef = reshape(elems(e, facePattern(:)), size(facePattern));
        n = size(ef, 1);
        faces(row:row+n-1, :) = ef;
        faceColor(row:row+n-1) = x(e);
        row = row + n;
    end

    patch('Vertices', model.mesh.nodes, 'Faces', faces, ...
        'FaceVertexCData', faceColor, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 1.0);
    xlabel('x'); ylabel('y'); zlabel('z');
end

function plotHistory(history, configs, resultRoot)
    fig = figure('Visible', 'off', 'Name', 'Test A convergence');
    tiledlayout(3, 1);

    nexttile;
    plot(history.iteration, history.J, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('J');

    nexttile;
    plot(history.iteration, history.C, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('Compliance C_k');
    legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), ...
        'Location', 'best');

    nexttile;
    plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('Volume fraction');

    saveas(fig, fullfile(resultRoot, 'convergence_history.png'));
    close(fig);
end

function plotLocationDensity(locationStats, resultRoot)
    loc = [locationStats.locationIndex]';
    meanDensity = [locationStats.meanDensity]';
    rhoGt05 = [locationStats.rhoGt05]';

    fig = figure('Visible', 'off', 'Name', 'Test A location density');
    tiledlayout(2, 1);

    nexttile;
    bar(loc, meanDensity);
    grid on; xlabel('Half-segment location'); ylabel('Mean density');

    nexttile;
    bar(loc, rhoGt05);
    grid on; xlabel('Half-segment location'); ylabel('Fraction \rho > 0.5');

    saveas(fig, fullfile(resultRoot, 'location_density.png'));
    close(fig);
end

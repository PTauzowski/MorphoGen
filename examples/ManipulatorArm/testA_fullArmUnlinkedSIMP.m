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
fdRelTol = 1.0e-2;   % p-norm aggregation + 48k elems: FD noise ~ 4-8e-3 expected for low-sensitivity elements

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
locationDensityRange = max([locationStats.meanDensity]) - min([locationStats.meanDensity]);
locationDensityStd = std([locationStats.meanDensity]);
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
    'history', 'configs', 'VolFrac', 'penal', 'pAgg', 'maxIter', ...
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

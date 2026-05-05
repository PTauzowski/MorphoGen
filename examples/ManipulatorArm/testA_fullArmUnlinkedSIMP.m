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
arm = armModelDefaults("thin");
E = arm.E;
nu = arm.nu;
R = arm.R;
r = arm.r;
h_seg = arm.h_seg;
alpha = arm.alpha;
res = arm.res;
res_th = arm.res_th;
Pz = arm.Pz;
ShapeFn = arm.ShapeFn;

VolFrac = 0.40;
penal = 3.0;
pAgg = 4.0;
maxIter = 500;
xminValue = arm.mma.xminValue;
Rfilter = arm.Rfilter;
useParallel = license('test', 'Distrib_Computing_Toolbox');

fdElemCount = 8;
fdStep = 1.0e-5;
fdRelTol = 1.0e-2;

moveLimit = arm.mma.moveLimit;
minMoveLimit = arm.mma.minMoveLimit;
moveDecay = arm.mma.moveDecay;
mmaDamping = arm.mma.mmaDamping;
changeTol = arm.mma.changeTol;
objectiveTol = 5.0e-3;
minIter = 10;

%configs = armLoadConfigs("sixPlusTension");
configs = armLoadConfigs("six");
nConfigs = numel(configs);
weights = ones(nConfigs, 1) / nConfigs;

resultRoot = fullfile(scriptDir, 'results', 'testA_fullArmUnlinkedSIMP');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test A full-arm unlinked SIMP diagnostic\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf('  VolFrac=%.3f, penal=%.2f, pAgg=%.2f, maxIter=%d, changeTol=%.1e, xmin=%.3f, parallel=%d\n', ...
    VolFrac, penal, pAgg, maxIter, changeTol, xminValue, useParallel);

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
        cfg.betas, ShapeFn, true, Pz, arm.constEndRing, arm.constMiddleRing);
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
const_elems = armConstRingElementIds(models{1}, arm, "full");
xmin(const_elems) = 1.0;
xmax(const_elems) = 1.0;
x(const_elems) = 1.0;
x = enforceVolumeFraction(x, VolFrac, xmin, xmax);
fprintf('\nConst ring elems: %d (end=%d, middle=%d)\n', ...
    numel(const_elems), arm.constEndRing, arm.constMiddleRing);

fprintf('\nBuilding full-arm sensitivity filter: Rfilter=%.6g, nDesign=%d\n', ...
    Rfilter, nDesign);
Wfilter = buildElementFilterMatrix(models{1}.mesh.nodes, models{1}.mesh.elems, ...
    (1:nDesign)', Rfilter);

%% ---- Compliance normalization and finite-difference gradient check -------
fprintf('\nComputing initial compliance normalization at x=VolFrac...\n');
[J0, grad0, C0, Cinit] = evaluateObjectiveAndGradient(analyses, x, penal, pAgg, weights, [], useParallel);
fprintf('  Initial J = %.8e\n', J0);
for k = 1:nConfigs
    fprintf('    %-12s C0 = %.8e\n', configs{k}.name, C0(k));
end

fprintf('\nFinite-difference gradient test (%d random elements, nominal h=%.1e)...\n', ...
    fdElemCount, fdStep);
fd = finiteDifferenceGradientTest(analyses, x, penal, pAgg, weights, C0, ...
    grad0, fdElemCount, fdStep, xmin, xmax, useParallel);
fprintf('  FD max relative error = %.3e\n', fd.maxRelativeError);
assert(fd.maxRelativeError < fdRelTol, ...
    'FD gradient check failed: max relative error %.3e exceeds %.3e.', ...
    fd.maxRelativeError, fdRelTol);

fdTable = struct2table(fd.rows);
writetable(fdTable, fullfile(resultRoot, 'finite_difference_gradient.csv'));

%% ---- MMA optimization ---------------------------------------------------
fprintf('\nRunning projected MMA for %d iterations...\n', maxIter);

opts = struct();
opts.penal = penal;
opts.pAgg = pAgg;
opts.weights = weights;
opts.C0 = C0;
opts.VolFrac = VolFrac;
opts.maxIter = maxIter;
opts.minIter = minIter;
opts.changeTol = changeTol;
opts.objectiveTol = objectiveTol;
opts.moveLimit = moveLimit;
opts.minMoveLimit = minMoveLimit;
opts.moveDecay = moveDecay;
opts.mmaDamping = mmaDamping;
opts.useParallel = useParallel;
opts.sensitivityFilter = Wfilter;
opts.fixedDesignVariables = const_elems;
opts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
opts.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

optResult = solveSIMPComplianceVolumeMMA(analyses, x, xmin, xmax, opts);

finalX = optResult.zFinal;
finalJ = optResult.finalJ;
finalC = optResult.finalC;
finalVolumeFraction = optResult.finalVolumeFraction;
finalChange = optResult.finalChange;
nIter = optResult.nIter;
history = optResult.history;
objectiveWindowConverged = optResult.objectiveWindowConverged;
objectiveDecreased = optResult.objectiveDecreased;
mmaAcceptable = optResult.mmaAcceptable;
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
    objectiveDecreased, history.J(1), finalJ);
fprintf('  Objective window converged   : %d\n', objectiveWindowConverged);
fprintf('  MMA acceptable               : %d\n', mmaAcceptable);
fprintf('  Volume constraint active     : %d (vf=%.6f)\n', volumeActive, finalVolumeFraction);
fprintf('  Topology nonuniform          : %d (std=%.4f, range=%.4f)\n', ...
    topologyNonuniform, std(finalX), max(finalX) - min(finalX));
fprintf('  Differs by arm location      : %d (half-seg mean range=%.4f)\n', ...
    topologyDiffersByArmLocation, locationDensityRange);

%% ---- Save histories and plots ------------------------------------------

save(fullfile(resultRoot, 'result.mat'), ...
    'finalX', 'finalJ', 'finalC', 'finalVolumeFraction', 'finalChange', ...
    'history', 'configs', 'nIter', 'VolFrac', 'penal', 'pAgg', 'maxIter', 'minIter', 'changeTol', ...
    'xminValue', 'Rfilter', 'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', ...
    'useParallel', 'const_elems', 'C0', 'Wfilter', ...
    'Pz', 'locationStats', 'locationDensityRange', 'locationDensityStd', ...
    'topologyDiffersByArmLocation', 'objectiveDecreased', ...
    'objectiveWindowConverged', 'mmaAcceptable', 'volumeActive', ...
    'topologyNonuniform');

saveHistoryCsv(resultRoot, history, configs);
saveSummaryCsv(resultRoot, fd, history, finalX, finalC, VolFrac, ...
    objectiveDecreased, objectiveWindowConverged, mmaAcceptable, volumeActive, ...
    topologyNonuniform, topologyDiffersByArmLocation, ...
    locationDensityRange, locationDensityStd);

%% ---- Post-processing: topology extraction ----------------------------------
fprintf('\nPost-processing topology extraction...\n');
ppOpts = struct();
ppOpts.penal        = penal;
ppOpts.pAgg         = pAgg;
ppOpts.weights      = weights;
ppOpts.VolFrac      = VolFrac;
ppOpts.fixedVars    = const_elems;
ppOpts.useParallel  = useParallel;
ppOpts.resultRoot   = resultRoot;
ppOpts.configNames  = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
%% ---- Structural performance metrics ----------------------------------------
fprintf('\nComputing structural performance metrics...\n');
x_ref_A = ones(nDesign, 1);
metrics_ref_A   = evaluateStructuralPerformance(analyses, x_ref_A, 1,     useParallel);
metrics_final_A = evaluateStructuralPerformance(analyses, finalX,  penal, useParallel);
saveStructuralMetricsCsv(metrics_ref_A, metrics_final_A, configs, finalVolumeFraction, resultRoot);

ppOpts.skipHeaviside = true;  % sensitivity-filtered SIMP: z IS physical density
ppOpts.runReanalysis = false;
ppOpts.runSweep     = false;

postResult = postprocessSIMPResult(optResult, models{1}, analyses, Wfilter, ppOpts);

plotFinalTopology(models{1}, finalX, resultRoot, postResult, analyses);
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

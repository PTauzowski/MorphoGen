% testH_fullArmStressLinkedSIMP
% Test H: full-arm linked stress-constrained SIMP formulation.

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(99, 'twister');

%% ---- Geometry and optimization parameters ----------------------------------
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

penal = 3.0;
maxIter = 200; 
xminValue = arm.mma.xminValue;
Rfilter = arm.Rfilter;
parallelWorkers = 6;
parallelAvailable = license('test', 'Distrib_Computing_Toolbox');
useParallel = parallelAvailable;
if useParallel
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('Processes', parallelWorkers);
    elseif pool.NumWorkers ~= parallelWorkers
        delete(pool);
        parpool('Processes', parallelWorkers);
    end
end

stressCoeff = 2.0;
stressPNorm = 16;
stressRelaxationQ = 0.5;
constraintTol = 5.0e-3;
activeConstraintTol = 2.0e-2;
changeTol = arm.mma.changeTol;
minIter = 10;
useProjection = true;
projectionBeta = 4.0;
projectionEta = 0.5;
nStressClusters = 3;

moveLimit = arm.mma.moveLimit;
minMoveLimit = 0.003;
moveDecay = 0.98;
mmaDamping = 0.5;

%configs = armLoadConfigs("sixPlusTension");
configs = armLoadConfigs("six");
nConfigs = numel(configs);

resultRoot = fullfile(scriptDir, 'results', 'testH_fullArmStressLinkedSIMP');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test H full-arm linked stress-constrained SIMP\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf(['  Objective=min volume, stressCoeff=%.3f, pStress=%.1f, qRelax=%.2f, ' ...
    'projection=%d(beta=%.1f, eta=%.2f), clusters=%d, maxIter=%d, parallel=%d, workers=%d\n'], ...
    stressCoeff, stressPNorm, stressRelaxationQ, useProjection, projectionBeta, projectionEta, ...
    nStressClusters, maxIter, useParallel, parallelWorkers);

%% ---- Build full-arm configurations -----------------------------------------
models = cell(nConfigs, 1);
analyses = cell(nConfigs, 1);
setupRows = cell(nConfigs, 1);

referenceElems = [];
referenceModel = [];
referenceElemCount = [];
referenceDofs = [];
referenceTaskDim = [];

for k = 1:nConfigs
    cfg = configs{k};
    fprintf('\nBuilding configuration %d/%d: %s, betas=%s\n', ...
        k, nConfigs, cfg.label, mat2str(cfg.betas));

    model = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
        cfg.betas, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
    analysis = model.analysis;
    nElems = analysis.getTotalElemsNumber();
    taskDim = analysis.getTaskDim();

    if k == 1
        referenceElems = model.mesh.elems;
        referenceModel = model;
        referenceElemCount = nElems;
        referenceDofs = analysis.ndofs;
        referenceTaskDim = taskDim;
        fprintf('  Reference element count : %d\n', referenceElemCount);
        fprintf('  Reference task DOFs     : %d\n', referenceTaskDim);
    else
        assert(nElems == referenceElemCount, 'Element-count mismatch in %s.', cfg.name);
        assert(taskDim == referenceTaskDim, 'DOF-count mismatch in %s.', cfg.name);
        assertLinkedArmLayoutCompatible(model, referenceModel, cfg.name);
        assert(isequal(analysis.ndofs, referenceDofs), 'DOF labels/order differ in %s.', cfg.name);
    end

    models{k} = model;
    analyses{k} = analysis;

    row.configName = string(cfg.name);
    row.configLabel = string(cfg.label);
    row.nNodes = size(model.mesh.nodes, 1);
    row.nElems = nElems;
    row.nTaskDofs = taskDim;
    row.nSupportedDofs = nnz(analysis.supports);
    row.nLoadedDofsBeforeSolve = nnz(analysis.Pnodal);
    row.sameConnectivityAsFirst = k == 1 || isequal(model.mesh.elems, referenceElems);
    row.sameLinkedLayoutAsFirst = true;
    row.sameDofsAsFirst = true;
    setupRows{k, 1} = row;
end

writetable(struct2table(vertcat(setupRows{:})), fullfile(resultRoot, 'configuration_setup.csv'));

%% ---- Linked design space ----------------------------------------------------
H = models{1}.halfSegmentNelems;
nElems = referenceElemCount;
nCopies = nElems / H;

fprintf('\nLinked design space: H=%d reference elements, nElems=%d, nCopies=%d\n', ...
    H, nElems, nCopies);
assert(mod(nCopies, 2) == 0, ...
    'nElems/H = %d is not even; segmentToArm requires an even number of copies.', nCopies);

xmin = xminValue * ones(H, 1);
xmax = ones(H, 1);
rho = ones(H, 1);
const_elems = armConstRingElementIds(models{1}, arm, "linked");
xmin(const_elems) = 1.0;
xmax(const_elems) = 1.0;
rho(const_elems) = 1.0;
fprintf('Const ring elems in linked rho: %d (end=%d, middle=%d)\n', ...
    numel(const_elems), arm.constEndRing, arm.constMiddleRing);

fprintf('\nBuilding linked reference sensitivity filter: Rfilter=%.6g, H=%d\n', ...
    Rfilter, H);
Wfilter = buildElementFilterMatrix(models{1}.mesh.nodes, models{1}.mesh.elems, ...
    (1:H)', Rfilter);

fprintf('\nBuilding and verifying linked sensitivity map...\n');
map = buildLinkedSensitivityMap(H, nElems);
e_ref_test = ceil(H / 2);
x_base = models{1}.segmentToArm(rho);
rhoPert = rho;
rhoPert(e_ref_test) = rhoPert(e_ref_test) + 1e-6;
x_pert = models{1}.segmentToArm(rhoPert);
changedElems = find(abs(x_pert - x_base) > 1e-10);
expectedAll = sort([map.normalMap(:, e_ref_test); map.flippedMap(:, e_ref_test)]);
assert(isequal(changedElems, expectedAll), ...
    'Linked map verification failed for e_ref=%d.', e_ref_test);
fprintf('  Map verification PASSED: e_ref=%d changes exactly %d elements\n', ...
    e_ref_test, numel(expectedAll));

%% ---- MMA optimization -------------------------------------------------------
fprintf('\nRunning stress-constrained MMA on rho (%d variables) for up to %d iterations...\n', ...
    H, maxIter);

opts = struct();
opts.penal = penal;
opts.stressCoeff = stressCoeff;
opts.stressPNorm = stressPNorm;
opts.stressRelaxationQ = stressRelaxationQ;
opts.maxIter = maxIter;
opts.minIter = minIter;
opts.changeTol = changeTol;
opts.constraintTol = constraintTol;
opts.activeConstraintTol = activeConstraintTol;
opts.moveLimit = moveLimit;
opts.minMoveLimit = minMoveLimit;
opts.moveDecay = moveDecay;
opts.mmaDamping = mmaDamping;
opts.useParallel = useParallel;
opts.sensitivityFilter = Wfilter;
opts.useProjection = useProjection;
opts.projectionBeta = projectionBeta;
opts.projectionEta = projectionEta;
opts.nStressClusters = nStressClusters;
opts.expand = @(z) models{1}.segmentToArm(z);
opts.pullback = @(dx) pullbackStressSet(dx, H, nElems);
opts.fixedDesignVariables = const_elems;
opts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
opts.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

optResult = solveSIMPVolumeStressMMA(analyses, rho, xmin, xmax, opts);

%% ---- Final/exported state and checks ---------------------------------------
rawRhoFinal = optResult.zFinal;
rawRhoPhysicalFinal = optResult.zPhysicalFinal;
rawXArmFinal = optResult.xFinal;
rawFinalStressAggregate = optResult.finalStressAggregate;
rawFinalMaxStress = optResult.finalMaxStress;
rawFinalVolumeFraction = optResult.finalVolumeFraction;
rawFinalConstraint = optResult.finalConstraint;
rawFinalChange = optResult.finalChange;
rawConstraintsSatisfied = optResult.constraintsSatisfied;
nIter = optResult.nIter;
history = optResult.history;
S0 = optResult.S0;
Starget = optResult.Starget;

exportIter = nIter;
usingBestFeasible = false;
if isfield(optResult, 'bestFeasible') && optResult.bestFeasible.valid
    exportIter = optResult.bestFeasible.iter;
    usingBestFeasible = exportIter ~= nIter;
end
exportIdx = exportIter + 1;

if usingBestFeasible
    rhoFinal = optResult.bestFeasible.z;
    rhoPhysicalFinal = optResult.bestFeasible.zPhysical;
    x_arm_final = opts.expand(rhoPhysicalFinal);
else
    rhoFinal = rawRhoFinal;
    rhoPhysicalFinal = rawRhoPhysicalFinal;
    x_arm_final = rawXArmFinal;
end

finalStressAggregate = history.stressAggregate(exportIdx, :)';
finalMaxStress = history.maxStress(exportIdx, :)';
finalVolumeFraction = mean(rhoPhysicalFinal);
finalConstraint = history.constraint(exportIdx, :)';
finalChange = history.change(exportIdx);
constraintsSatisfied = max(finalConstraint) <= constraintTol;
rhoNonuniform = std(rhoPhysicalFinal) > 0.03 && (max(rhoPhysicalFinal) - min(rhoPhysicalFinal)) > 0.15;
topologyNonuniform = rhoNonuniform;
locationDensityRange = max(rhoPhysicalFinal) - min(rhoPhysicalFinal);
locationDensityStd = std(rhoPhysicalFinal);
topologyDiffersByArmLocation = locationDensityRange > 0.03;

locationStats = computeLocationDensityStats(models{1}, x_arm_final);
writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

fprintf('\nFinal checks\n');
fprintf('  Exported iteration          : %d / %d (best feasible=%d)\n', ...
    exportIter, nIter, usingBestFeasible);
fprintf('  Final volume fraction        : %.6f\n', finalVolumeFraction);
fprintf('  Stress constraints ok        : %d (max g=%.3e)\n', constraintsSatisfied, max(finalConstraint));
fprintf('  Rho nonuniform               : %d (std=%.4f, range=%.4f)\n', ...
    rhoNonuniform, std(rhoPhysicalFinal), max(rhoPhysicalFinal) - min(rhoPhysicalFinal));

%% ---- Save histories and plots ----------------------------------------------
save(fullfile(resultRoot, 'result.mat'), ...
    'rhoFinal', 'rhoPhysicalFinal', 'x_arm_final', 'finalStressAggregate', 'finalMaxStress', ...
    'finalVolumeFraction', 'finalConstraint', 'finalChange', 'history', ...
    'configs', 'map', 'H', 'nElems', 'nCopies', 'nIter', ...
    'exportIter', 'usingBestFeasible', 'optResult', ...
    'rawRhoFinal', 'rawRhoPhysicalFinal', 'rawXArmFinal', ...
    'rawFinalStressAggregate', 'rawFinalMaxStress', 'rawFinalVolumeFraction', ...
    'rawFinalConstraint', 'rawFinalChange', 'rawConstraintsSatisfied', ...
    'S0', 'Starget', 'stressCoeff', 'stressPNorm', 'stressRelaxationQ', ...
    'useProjection', 'projectionBeta', 'projectionEta', 'nStressClusters', ...
    'penal', 'maxIter', 'minIter', 'changeTol', 'constraintTol', 'activeConstraintTol', ...
    'xminValue', 'Rfilter', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', ...
    'useParallel', 'parallelAvailable', 'parallelWorkers', 'Pz', ...
    'const_elems', 'locationStats', 'locationDensityRange', 'locationDensityStd', ...
    'topologyDiffersByArmLocation', 'constraintsSatisfied', 'topologyNonuniform', 'rhoNonuniform');

saveStressSIMPHistoryCsv(resultRoot, history, configs);
saveStressSIMPSummaryCsv(resultRoot, history, rhoPhysicalFinal, finalStressAggregate, ...
    finalMaxStress, finalConstraint, stressCoeff, stressPNorm, stressRelaxationQ, ...
    constraintsSatisfied, topologyNonuniform, topologyDiffersByArmLocation, ...
    locationDensityRange, locationDensityStd);

%% ---- Post-processing: topology extraction ----------------------------------
fprintf('\nPost-processing topology extraction...\n');
ppOpts = struct();
ppOpts.penal              = penal;
ppOpts.stressPNorm        = stressPNorm;
ppOpts.stressRelaxationQ  = stressRelaxationQ;
ppOpts.nStressClusters    = nStressClusters;
ppOpts.Starget            = Starget;
ppOpts.constraintTol      = constraintTol;
ppOpts.fixedVars          = const_elems;
ppOpts.useParallel        = useParallel;
ppOpts.resultRoot         = resultRoot;
ppOpts.configNames        = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
ppOpts.plotModels         = models;
ppOpts.plotConfigs        = configs;
ppOpts.runReanalysis      = false;

ppOptResult = optResult;
ppOptResult.zFinal = rhoFinal;
ppOptResult.zPhysicalFinal = rhoPhysicalFinal;
ppOptResult.xFinal = x_arm_final;
ppOptResult.finalVolumeFraction = finalVolumeFraction;
ppOptResult.finalStressAggregate = finalStressAggregate;
ppOptResult.finalMaxStress = finalMaxStress;
ppOptResult.finalConstraint = finalConstraint;
postResult = postprocessStressSIMPResult(ppOptResult, models{1}, analyses, Wfilter, ppOpts);

plotFinalTopology(models{1}, x_arm_final, resultRoot, postResult, analyses);
plotTopologyConfigurations(models, x_arm_final > 0.5, resultRoot, ...
    "final_threshold_rho_gt_05_by_config", "Final rho > 0.5 by configuration", configs);
plotStressSIMPHistory(history, configs, resultRoot, 'Test H');
plotLocationDensity(locationStats, resultRoot);

%% ---- Structural performance metrics ----------------------------------------
fprintf('\nComputing structural performance metrics...\n');
x_ref_H = ones(nElems, 1);
metrics_ref_H   = evaluateStructuralPerformance(analyses, x_ref_H,    1,     useParallel);
metrics_final_H = evaluateStructuralPerformance(analyses, x_arm_final, penal, useParallel);
saveStructuralMetricsCsv(metrics_ref_H, metrics_final_H, configs, finalVolumeFraction, resultRoot);

fig = figure('Name', 'Test H reference rho');
bar(rhoPhysicalFinal, 'FaceColor', [0.28 0.45 0.72]);
xlabel('Reference half-segment element index');
ylabel('Density \rho');
title(sprintf('Test H: reference half-segment density (H=%d, vf=%.3f)', ...
    H, mean(rhoPhysicalFinal)));
ylim([0 1.05]);
exportgraphics(fig, fullfile(resultRoot, 'rho_final.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'rho_final.fig'));
close(fig);

fprintf('\nSaved Test H outputs to %s\n', resultRoot);

if ~constraintsSatisfied
    warning('TestH:StressConstraintsViolated', ...
        'Final stress constraints are violated: max g=%.3e.', max(finalConstraint));
end

%% ---- Local helpers ----------------------------------------------------------
function dSdrho = pullbackStressSet(dSdx, H, nElems)
    nConfigs = size(dSdx, 2);
    dSdrho = zeros(H, nConfigs);
    for k = 1:nConfigs
        dSdrho(:, k) = pullbackFullArmSensitivity(dSdx(:, k), H, nElems);
    end
end

% testG_fullArmStressUnlinkedSIMP
% Test G: full-arm unlinked stress-constrained SIMP formulation.

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(88, 'twister');

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

%% ---- Robust beta options (set robustEnabled=true to activate CG robust loop) --
robustEnabled     = false;
robustDeltaDeg    = 90;
robustStressRatio = 3.0;    % stressLimit = ratio × full-pipe max-HM stress
robustDispRatio   = Inf;    % Inf = displacement constraint inactive
robustMaxCGIter   = 5;
robustTopK        = 3;
robustPropLevel   = 1;

if robustEnabled
    resultRoot = fullfile(scriptDir, 'results', ...
        sprintf('testG_fullArmStressUnlinkedSIMP_robust_d%g', robustDeltaDeg));
else
    resultRoot = fullfile(scriptDir, 'results', 'testG_fullArmStressUnlinkedSIMP');
end
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test G full-arm unlinked stress-constrained SIMP\n');
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
        referenceElemCount = nElems;
        referenceDofs = analysis.ndofs;
        referenceTaskDim = taskDim;
        fprintf('  Reference element count: %d\n', referenceElemCount);
        fprintf('  Reference task DOFs    : %d\n', referenceTaskDim);
    else
        assert(nElems == referenceElemCount, 'Element-count mismatch in %s.', cfg.name);
        assert(taskDim == referenceTaskDim, 'DOF-count mismatch in %s.', cfg.name);
        assert(isequal(model.mesh.elems, referenceElems), 'Mesh connectivity differs in %s.', cfg.name);
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
    setupRows{k, 1} = row;
end

writetable(struct2table(vertcat(setupRows{:})), fullfile(resultRoot, 'configuration_setup.csv'));

%% ---- Design space -----------------------------------------------------------
nDesign = referenceElemCount;
xmin = xminValue * ones(nDesign, 1);
xmax = ones(nDesign, 1);
x = ones(nDesign, 1);
const_elems = armConstRingElementIds(models{1}, arm, "full");
xmin(const_elems) = 1.0;
xmax(const_elems) = 1.0;
x(const_elems) = 1.0;
fprintf('\nConst ring elems: %d (end=%d, middle=%d)\n', ...
    numel(const_elems), arm.constEndRing, arm.constMiddleRing);

fprintf('\nBuilding full-arm sensitivity filter: Rfilter=%.6g, nDesign=%d\n', ...
    Rfilter, nDesign);
Wfilter = buildElementFilterMatrix(models{1}.mesh.nodes, models{1}.mesh.elems, ...
    (1:nDesign)', Rfilter);

%% ---- MMA optimization -------------------------------------------------------
fprintf('\nRunning stress-constrained MMA on x (%d variables) for up to %d iterations...\n', ...
    nDesign, maxIter);

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
opts.fixedDesignVariables = const_elems;
opts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
opts.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

if robustEnabled
    fprintf('\nComputing full-pipe reference for robust stress/disp limits...\n');
    evalOpts.penal = 1;
    [fpMetrics, ~] = evaluateLinkedDensityMetrics(analyses, ones(nDesign, 1), evalOpts);
    stressLimitRobust = robustStressRatio * max(fpMetrics.maxHM);
    dispLimitRobust   = robustDispRatio   * max(abs(fpMetrics.tipUz));
    fprintf('  Full-pipe: maxHM=%.4e Pa, |tipUz|=%.4e m\n', ...
        max(fpMetrics.maxHM), max(abs(fpMetrics.tipUz)));
    fprintf('  Stress limit: %.4e Pa (x%.2f)  Disp limit: %.4e m (x%.2f)\n', ...
        stressLimitRobust, robustStressRatio, dispLimitRobust, robustDispRatio);

    cgOpts.maxCGIter  = robustMaxCGIter;
    cgOpts.topK       = robustTopK;
    cgOpts.propLevel  = robustPropLevel;
    cgOpts.penal      = penal;
    cgOpts.verbose    = true;

    [optResult, adversarialConfigs, cgHistory] = constraintGenerationRobustSIMP( ...
        analyses, models{1}, arm, opts, stressLimitRobust, dispLimitRobust, ...
        x, xmin, xmax, robustDeltaDeg, cgOpts);
else
    optResult          = solveSIMPVolumeStressMMA(analyses, x, xmin, xmax, opts);
    adversarialConfigs = struct([]);
    cgHistory          = struct([]);
end

%% ---- Final/exported state and checks ---------------------------------------
rawFinalX = optResult.zFinal;
rawFinalXPhysical = optResult.zPhysicalFinal;
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
    finalX = optResult.bestFeasible.z;
    finalXPhysical = optResult.bestFeasible.zPhysical;
else
    finalX = rawFinalX;
    finalXPhysical = rawFinalXPhysical;
end

finalStressAggregate = history.stressAggregate(exportIdx, :)';
finalMaxStress = history.maxStress(exportIdx, :)';
finalVolumeFraction = mean(finalXPhysical);
finalConstraint = history.constraint(exportIdx, :)';
finalChange = history.change(exportIdx);
constraintsSatisfied = max(finalConstraint) <= constraintTol;
topologyNonuniform = std(finalXPhysical) > 0.03 && (max(finalXPhysical) - min(finalXPhysical)) > 0.15;

locationStats = computeLocationDensityStats(models{1}, finalXPhysical);
locationDensityRange = max([locationStats.meanDensity]) - min([locationStats.meanDensity]);
locationDensityStd = std([locationStats.meanDensity]);
topologyDiffersByArmLocation = locationDensityRange > 0.03;
writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

fprintf('\nFinal checks\n');
fprintf('  Exported iteration          : %d / %d (best feasible=%d)\n', ...
    exportIter, nIter, usingBestFeasible);
fprintf('  Final volume fraction        : %.6f\n', finalVolumeFraction);
fprintf('  Stress constraints ok        : %d (max g=%.3e)\n', constraintsSatisfied, max(finalConstraint));
fprintf('  Topology nonuniform          : %d (std=%.4f, range=%.4f)\n', ...
    topologyNonuniform, std(finalXPhysical), max(finalXPhysical) - min(finalXPhysical));

%% ---- Save histories and plots ----------------------------------------------
save(fullfile(resultRoot, 'result.mat'), ...
    'finalX', 'finalXPhysical', 'finalStressAggregate', 'finalMaxStress', 'finalVolumeFraction', ...
    'finalConstraint', 'finalChange', 'history', 'configs', 'nIter', ...
    'exportIter', 'usingBestFeasible', 'optResult', ...
    'rawFinalX', 'rawFinalXPhysical', 'rawFinalStressAggregate', 'rawFinalMaxStress', ...
    'rawFinalVolumeFraction', 'rawFinalConstraint', 'rawFinalChange', 'rawConstraintsSatisfied', ...
    'S0', 'Starget', 'stressCoeff', 'stressPNorm', 'stressRelaxationQ', ...
    'useProjection', 'projectionBeta', 'projectionEta', 'nStressClusters', ...
    'penal', 'maxIter', 'minIter', 'changeTol', 'constraintTol', 'activeConstraintTol', ...
    'xminValue', 'Rfilter', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', ...
    'useParallel', 'parallelAvailable', 'parallelWorkers', 'Pz', ...
    'const_elems', 'locationStats', 'locationDensityRange', 'locationDensityStd', ...
    'topologyDiffersByArmLocation', 'constraintsSatisfied', 'topologyNonuniform');

saveStressSIMPHistoryCsv(resultRoot, history, configs);
saveStressSIMPSummaryCsv(resultRoot, history, finalXPhysical, finalStressAggregate, ...
    finalMaxStress, finalConstraint, stressCoeff, stressPNorm, stressRelaxationQ, ...
    constraintsSatisfied, topologyNonuniform, topologyDiffersByArmLocation, ...
    locationDensityRange, locationDensityStd);

plotFinalTopology(models{1}, finalXPhysical, resultRoot, [], analyses);
plotStressSIMPHistory(history, configs, resultRoot, 'Test G');
plotLocationDensity(locationStats, resultRoot);

%% ---- Structural performance metrics ----------------------------------------
fprintf('\nComputing structural performance metrics...\n');
x_ref_G = ones(nDesign, 1);
metrics_ref_G   = evaluateStructuralPerformance(analyses, x_ref_G,       1,     useParallel);
metrics_final_G = evaluateStructuralPerformance(analyses, finalXPhysical, penal, useParallel);
saveStructuralMetricsCsv(metrics_ref_G, metrics_final_G, configs, finalVolumeFraction, resultRoot);

if robustEnabled && ~isempty(adversarialConfigs)
    save(fullfile(resultRoot, 'robust_result.mat'), ...
        'adversarialConfigs', 'cgHistory', 'stressLimitRobust', 'dispLimitRobust', ...
        'robustDeltaDeg', 'robustStressRatio', 'robustDispRatio');
    rows = struct([]);
    for ci = 1:numel(adversarialConfigs)
        c = adversarialConfigs(ci);
        row.ci = ci;  row.stressRatio = c.stressRatio;  row.dispRatio = c.dispRatio;
        row.violated = double(c.violated);
        row.maxHM_MPa = c.maxHM / 1e6;  row.tipUz_m = c.tipUz;
        bv = c.beta;
        for ji = 1:numel(bv), row.(sprintf('beta%d', ji)) = bv(ji); end
        rows = [rows; row]; %#ok<AGROW>
    end
    if ~isempty(rows)
        writetable(struct2table(rows), fullfile(resultRoot, 'adversarial_config_summary.csv'));
    end
end

fprintf('\nSaved Test G outputs to %s\n', resultRoot);

if ~constraintsSatisfied
    warning('TestG:StressConstraintsViolated', ...
        'Final stress constraints are violated: max g=%.3e.', max(finalConstraint));
end

% testD_fullArmUnlinkedStressIntensity
% Test D: full-arm unlinked multi-configuration stress-intensity ESO.
%
% Formulation:
%   Design variables: x [nElems x 1] -- one independent density per full-arm
%   solid element.  This is the stress-intensity counterpart of Test A and
%   does not enforce repeated-module linking.
%
%   Set stressAggregation="max" for the max-over-configs stress-intensity
%   envelope or "average" for summed/average intensity.

clear; close all; clc;
delete(gcp('nocreate'));
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

parpool('Processes', 6);

rng(44, 'twister');

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

VolFrac   = 0.40;
penal     = 3.0;
maxIter   = 200;
xminValue = arm.mma.xminValue;

useParallel = true;
parallelAvailable = license('test', 'Distrib_Computing_Toolbox');

% Optional second-order / buckling-amplified analysis.
useSecondOrderBuckling = false;
bucklingLambdaFactor = 0.95;

% Stress-intensity aggregation over load configurations: "max" or "average".
stressAggregation = "max";

% Stress-intensity ESO parameters
Rmin    = arm.Rfilter;
maxais  = 0.03;
max_elem_removal_factor = 0.05;

%configs = armLoadConfigs("sixPlusTension");
configs = armLoadConfigs("six");
nConfigs = numel(configs);

resultRoot = fullfile(scriptDir, 'results', 'testD_fullArmUnlinkedStressIntensity');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test D full-arm unlinked stress-intensity ESO\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf('  VolFrac=%.3f, penal=%.2f, maxIter=%d, Rmin=%.4f, maxais=%.4f, parallel=%d, parallelAvailable=%d\n', ...
    VolFrac, penal, maxIter, Rmin, maxais, useParallel, parallelAvailable);
fprintf('  Second-order buckling effects: %d (lambda factor=%.3f)\n', ...
    useSecondOrderBuckling, bucklingLambdaFactor);
fprintf('  Stress-intensity aggregation: %s\n', stressAggregation);

%% ---- Build full-arm configurations -----------------------------------------
models    = cell(nConfigs, 1);
analyses  = cell(nConfigs, 1);
setupRows = cell(nConfigs, 1);

referenceElems     = [];
referenceElemCount = [];
referenceDofs      = [];
referenceTaskDim   = [];

for k = 1:nConfigs
    cfg = configs{k};
    fprintf('\nBuilding configuration %d/%d: %s, betas=%s\n', ...
        k, nConfigs, cfg.label, mat2str(cfg.betas));

    model = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
        cfg.betas, ShapeFn, true, Pz, arm.constEndRing, arm.constMiddleRing);
    analysis = model.analysis;
    if useSecondOrderBuckling
        analysis = SecondOrderElasticityWeighted(model.fe, model.mesh, bucklingLambdaFactor, false);
        analysis.Pnodal = model.analysis.Pnodal;
        analysis.Pfem = model.analysis.Pfem;
        analysis.supports = model.analysis.supports;
        analysis.rotations = model.analysis.rotations;
    end

    nElems = analysis.getTotalElemsNumber();
    taskDim = analysis.getTaskDim();
    nSupports = nnz(analysis.supports);
    nLoadedDofs = nnz(analysis.Pnodal);

    if k == 1
        referenceElems = model.mesh.elems;
        referenceElemCount = nElems;
        referenceDofs = analysis.ndofs;
        referenceTaskDim = taskDim;
        fprintf('  Reference element count : %d\n', referenceElemCount);
        fprintf('  Reference task DOFs     : %d\n', referenceTaskDim);
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
    row.analysisClass = string(class(analysis));
    row.useSecondOrderBuckling = useSecondOrderBuckling;
    row.bucklingLambdaFactor = bucklingLambdaFactor;
    setupRows{k, 1} = row;

    fprintf('  nodes=%d, elems=%d, taskDOFs=%d, supportedDOFs=%d, loadedDOFs=%d\n', ...
        row.nNodes, row.nElems, row.nTaskDofs, row.nSupportedDofs, row.nLoadedDofsBeforeSolve);
end

setupTable = struct2table(vertcat(setupRows{:}));
writetable(setupTable, fullfile(resultRoot, 'configuration_setup.csv'));

nElems = referenceElemCount;
fprintf('\nUnlinked design space: nElems=%d full-arm design variables\n', nElems);
const_elems = armConstRingElementIds(models{1}, arm, "full");
fprintf('Const ring elems: %d (end=%d, middle=%d)\n', ...
    numel(const_elems), arm.constEndRing, arm.constMiddleRing);

%% ---- Create stress-intensity topology optimization object -------------------
fprintf('\nBuilding stress-intensity topology optimization (filter matrix)...\n');
analysesArray = [analyses{:}];
switch lower(stressAggregation)
    case "max"
        topOpt = StressIntensityMultiMaxTopologyOptimization(Rmin, analysesArray, maxais, penal, VolFrac, false);
    case {"average", "avg", "mean"}
        topOpt = StressIntensityMultiAvTopologyOptimization(Rmin, analysesArray, maxais, penal, VolFrac, false);
        stressAggregation = "average";
    otherwise
        error('Unknown stressAggregation "%s". Use "max" or "average".', stressAggregation);
end
topOpt.allx = [];
topOpt.useParallel = useParallel && parallelAvailable;
topOpt.setConstElems(const_elems);
fprintf('  Done. Filter weights: [%d x %d] sparse.\n', ...
    size(topOpt.weights, 1), size(topOpt.weights, 2));
fprintf('  Aggregation class: %s\n', class(topOpt));

%% ---- Custom unlinked ESO loop -----------------------------------------------
fprintf('\nRunning unlinked ESO for up to %d iterations (target vf=%.3f)...\n', maxIter, VolFrac);

x = ones(nElems, 1);
erased_x = false(nElems, 1);
x(const_elems) = 1.0;

xHistory = zeros(nElems, maxIter);
volHistory = nan(maxIter, 1);
erasedHistory = nan(maxIter, 1);
maxStressHistory = nan(maxIter, nConfigs);
maxDisplacementHistory = nan(maxIter, nConfigs);
iterTimeHistory = nan(maxIter, 1);

nIter = 0;
for iter = 1:maxIter
    iterTic = tic;

    topOpt.x = x;

    nStressBefore = numel(topOpt.maxstress);
    nDisplacementBefore = numel(topOpt.maxdisplacement);
    ais_full = topOpt.weights * topOpt.computeAverageIntensities();
    maxStressHistory(iter, :) = topOpt.maxstress(nStressBefore + 1 : end);
    maxDisplacementHistory(iter, :) = topOpt.maxdisplacement(nDisplacementBefore + 1 : end);

    notErasedID = setdiff(find(~erased_x), const_elems);
    if ~isempty(notErasedID)
        ais_ne = ais_full(notErasedID);
        ais_range = max(ais_ne) - min(ais_ne);
        removeList = ais_ne < min(ais_ne) + ais_range * maxais;
        maxRemove = round(nElems * max_elem_removal_factor);
        if sum(removeList) > maxRemove
            [~, ai] = sort(ais_ne);
            removeList = false(numel(notErasedID), 1);
            removeList(ai(1:maxRemove)) = true;
        end
        erased_x(notErasedID(removeList)) = true;
    end

    x(erased_x) = min(1, max(xminValue, x(erased_x) .* ais_full(erased_x) .^ penal));
    x(~erased_x) = 1;
    x(const_elems) = 1.0;

    vol = mean(x);
    iterTimeSec = toc(iterTic);

    xHistory(:, iter) = x;
    volHistory(iter) = vol;
    erasedHistory(iter) = sum(erased_x);
    iterTimeHistory(iter) = iterTimeSec;
    nIter = iter;

    fprintf('%4d  vol=%.4f  erased=%d/%d  time=%.2fs  maxStress=[', ...
        iter, vol, sum(erased_x), nElems, iterTimeSec);
    fprintf('%.3e ', maxStressHistory(iter, :));
    fprintf(']\n');

    if vol <= VolFrac
        fprintf('  Target volume fraction %.3f reached at iteration %d.\n', VolFrac, iter);
        break;
    end
end

%% ---- Final state and pass checks -------------------------------------------
xFinal = x;
finalVolumeFraction = mean(xFinal);
nErasedFinal = sum(erased_x);

volumeReached = finalVolumeFraction <= VolFrac + 1e-2;
xNonuniform = std(xFinal) > 0.03 && (max(xFinal) - min(xFinal)) > 0.15;
someErased = nErasedFinal > 0;

fprintf('\nFinal checks\n');
fprintf('  Volume fraction reached   : %d (vf=%.6f, target=%.3f)\n', ...
    volumeReached, finalVolumeFraction, VolFrac);
fprintf('  X nonuniform              : %d (std=%.4f, range=%.4f)\n', ...
    xNonuniform, std(xFinal), max(xFinal) - min(xFinal));
fprintf('  Elements erased in x      : %d / %d (%.1f%%)\n', ...
    nErasedFinal, nElems, 100 * nErasedFinal / nElems);

%% ---- Location-density diagnostics on full-arm topology ---------------------
locationStats = computeLocationDensityStats(models{1}, xFinal);
writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

%% ---- Save results ----------------------------------------------------------
xHistory = xHistory(:, 1:nIter);
volHistory = volHistory(1:nIter);
erasedHistory = erasedHistory(1:nIter);
maxStressHistory = maxStressHistory(1:nIter, :);
maxDisplacementHistory = maxDisplacementHistory(1:nIter, :);
iterTimeHistory = iterTimeHistory(1:nIter);

history.iteration = (1:nIter)';
history.volumeFraction = volHistory;
history.nErased = erasedHistory;
history.iterationTimeSec = iterTimeHistory;
history.maxStress = maxStressHistory;
history.maxDisplacement = maxDisplacementHistory;
history.x = xHistory;
history.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
history.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

save(fullfile(resultRoot, 'result.mat'), ...
    'xFinal', 'finalVolumeFraction', 'nErasedFinal', ...
    'erased_x', 'history', 'configs', 'const_elems', ...
    'nElems', 'nIter', ...
    'VolFrac', 'penal', 'Rmin', 'maxais', 'maxIter', 'xminValue', ...
    'useParallel', 'parallelAvailable', 'useSecondOrderBuckling', 'bucklingLambdaFactor', ...
    'stressAggregation', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', 'Pz', ...
    'locationStats', 'xNonuniform', 'volumeReached', 'someErased');

histRows = struct();
for i = 1:nIter
    histRows(i).iteration = i;
    histRows(i).volumeFraction = volHistory(i);
    histRows(i).nErased = erasedHistory(i);
    histRows(i).iterationTimeSec = iterTimeHistory(i);
    for k = 1:nConfigs
        histRows(i).(['maxStress_' configs{k}.name]) = maxStressHistory(i, k);
        histRows(i).(['maxDisplacement_' configs{k}.name]) = maxDisplacementHistory(i, k);
    end
end
writetable(struct2table(histRows), fullfile(resultRoot, 'history.csv'));

%% ---- Post-processing: optimal iteration selection --------------------------
fprintf('\nPost-processing: selecting optimal ESO iteration...\n');
ppOpts = struct();
ppOpts.VolFrac       = VolFrac;
% expand is identity for unlinked: history.x is already full-arm
ppOpts.penal         = penal;
ppOpts.useParallel   = useParallel;
ppOpts.resultRoot    = resultRoot;
ppOpts.configNames   = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
ppOpts.runReanalysis = false;

postResult = postprocessStressIntensityResult(history, models{1}, analyses, ppOpts);

plotFinalTopology(models{1}, xFinal, resultRoot, postResult, analyses);
plotLocationDensity(locationStats, resultRoot);

%% ---- Structural performance metrics ----------------------------------------
fprintf('\nComputing structural performance metrics...\n');
x_ref_D = ones(nElems, 1);
metrics_ref_D   = evaluateStructuralPerformance(analyses, x_ref_D, 1,     useParallel);
metrics_final_D = evaluateStructuralPerformance(analyses, xFinal,  penal, useParallel);
saveStructuralMetricsCsv(metrics_ref_D, metrics_final_D, configs, finalVolumeFraction, resultRoot);

fig = figure('Name', 'Test D volume history');
plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.2);
yline(VolFrac, '--r', sprintf('Target %.2f', VolFrac));
grid on;
xlabel('Iteration'); ylabel('Volume fraction');
title('Test D: unlinked ESO volume fraction');
saveas(fig, fullfile(resultRoot, 'volume_history.png'));
savefig(fig, fullfile(resultRoot, 'volume_history.fig'));
close(fig);

fig = figure('Name', 'Test D max stress history');
plot(history.iteration, history.maxStress, '-o', 'LineWidth', 1.0);
grid on;
xlabel('Iteration'); ylabel('Max HM stress');
legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), 'Location', 'best');
title('Test D: per-config max stress vs iteration');
saveas(fig, fullfile(resultRoot, 'stress_history.png'));
savefig(fig, fullfile(resultRoot, 'stress_history.fig'));
close(fig);

fig = figure('Name', 'Test D max displacement history');
plot(history.iteration, history.maxDisplacement, '-o', 'LineWidth', 1.0);
grid on;
xlabel('Iteration'); ylabel('Max displacement magnitude');
legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), 'Location', 'best');
title('Test D: per-config max displacement vs iteration');
saveas(fig, fullfile(resultRoot, 'displacement_history.png'));
savefig(fig, fullfile(resultRoot, 'displacement_history.fig'));
close(fig);

fig = figure('Name', 'Test D full-arm x');
bar(xFinal, 'FaceColor', [0.42 0.42 0.42]);
xlabel('Full-arm element index');
ylabel('Density x');
title(sprintf('Test D: full-arm x (nElems=%d, vf=%.3f)', nElems, mean(xFinal)));
ylim([0 1.05]);
exportgraphics(fig, fullfile(resultRoot, 'x_final.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'x_final.fig'));
close(fig);

%% ---- Optional: compare to Test C (linked ESO) and Test A (unlinked MMA) -----
testCPath = fullfile(scriptDir, 'results', 'testC_fullArmLinkedStressIntensity', 'result.mat');
testAPath = fullfile(scriptDir, 'results', 'testA_fullArmUnlinkedSIMP', 'result.mat');
if exist(testCPath, 'file') || exist(testAPath, 'file')
    fprintf('\nComparison\n');
    compRow.testDVolFrac = finalVolumeFraction;
    compRow.testDnIter = nIter;
    compRow.testDnErased = nErasedFinal;
    if exist(testCPath, 'file')
        resC = load(testCPath, 'finalVolumeFraction', 'nIter', 'nErasedFinal');
        fprintf('  Test C (linked ESO)    vf=%.6f  nIter=%d  erased=%d\n', ...
            resC.finalVolumeFraction, resC.nIter, resC.nErasedFinal);
        compRow.testCVolFrac = resC.finalVolumeFraction;
        compRow.testCnIter = resC.nIter;
        compRow.testCnErased = resC.nErasedFinal;
    end
    if exist(testAPath, 'file')
        resA = load(testAPath, 'finalJ', 'finalVolumeFraction');
        fprintf('  Test A (unlinked MMA)  vf=%.6f  J=%.6e\n', ...
            resA.finalVolumeFraction, resA.finalJ);
        compRow.testAFinalJ = resA.finalJ;
        compRow.testAVolFrac = resA.finalVolumeFraction;
    end
    writetable(struct2table(compRow), fullfile(resultRoot, 'comparison.csv'));
end

fprintf('\nSaved Test D outputs to %s\n', resultRoot);

%% ---- Pass criteria warnings ------------------------------------------------
if ~volumeReached
    warning('TestD:VolumeTargetNotReached', ...
        'Volume fraction %.6f did not reach target %.3f.', finalVolumeFraction, VolFrac);
end
if ~xNonuniform
    warning('TestD:XTooUniform', ...
        'Full-arm x is too uniform after ESO: std=%.4f, range=%.4f.', ...
        std(xFinal), max(xFinal) - min(xFinal));
end
if ~someErased
    warning('TestD:NoElementsErased', ...
        'No x elements were erased -- ESO loop did not run or remove anything.');
end

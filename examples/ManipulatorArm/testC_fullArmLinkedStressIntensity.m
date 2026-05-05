% testC_fullArmLinkedStressIntensity
% Test C: full-arm linked multi-configuration stress-intensity envelope (ESO).
%
% Formulation:
%   Design variables: rho [H x 1] — reference half-segment densities.
%   Full-arm expansion: x_arm = model.segmentToArm(rho)
%   Pattern: [rho; flip(rho); rho; flip(rho); ...]
%
%   The optimizer uses a stress-intensity approach instead of the
%   compliance/MMA approach in testB.  Set stressAggregation="max" for the
%   max-over-configs envelope or "average" for summed/average intensity.
%
%   Each iteration:
%     1.  Expand rho -> x_arm.
%     2.  Solve all solid FE problems with x_arm.^penal stiffness.
%     3.  Compute element von Mises (HM) stress, take max over configs.
%     4.  Apply spatial filter -> ais_full [nElems x 1].
%     5.  Pull back to rho-space: ais_rho = mean over copies.
%     6.  Mark the least-stressed rho elements as erased (ESO rule).
%     7.  Scale erased rho entries by (ais_rho / max)^penal;
%         un-erased entries stay at 1.
%     8.  Stop when mean(rho) <= VolFrac.
%
% Configurations (3):
%   1. Max M_z bending
%   2. Max M_s torsion
%   3. Max T_y shear
%
% Pass criteria:
%   - Volume fraction reached target.
%   - Reference module rho is nonuniform.
%   - Some elements erased (sum(erased_rho) > 0).

clear; close all; clc;
delete(gcp('nocreate'));
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

parpool('Processes', 6);

rng(33, 'twister');

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

useParallel = true;  % Set true to use parfor over load configurations.
parallelAvailable = license('test', 'Distrib_Computing_Toolbox');

% Optional second-order / buckling-amplified analysis, following the
% SecondOrderElasticityWeighted use in configurationTopologyMulti.
useSecondOrderBuckling = false;
bucklingLambdaFactor = 0.95;

% Stress-intensity aggregation over load configurations: "max" or "average".
stressAggregation = "max";

% Stress-intensity ESO parameters
Rmin    = arm.Rfilter;
maxais  = 0.03;        % fraction of ais range to remove per step
max_elem_removal_factor = 0.05;  % hard cap: max fraction of H removed per step

%configs = armLoadConfigs("sixPlusTension");
configs = armLoadConfigs("six");

nConfigs = numel(configs);

resultRoot = fullfile(scriptDir, 'results', 'testC_fullArmLinkedStressIntensity');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test C full-arm linked stress-intensity ESO\n');
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
referenceModel     = [];
referenceElemCount = [];
referenceDofs      = [];
referenceTaskDim   = [];

for k = 1:nConfigs
    cfg = configs{k};
    fprintf('\nBuilding configuration %d/%d: %s, betas=%s\n', ...
        k, nConfigs, cfg.label, mat2str(cfg.betas));

    model    = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
        cfg.betas, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
    analysis = model.analysis;
    if useSecondOrderBuckling
        analysis = SecondOrderElasticityWeighted(model.fe, model.mesh, bucklingLambdaFactor, false);
        analysis.Pnodal = model.analysis.Pnodal;
        analysis.Pfem = model.analysis.Pfem;
        analysis.supports = model.analysis.supports;
        analysis.rotations = model.analysis.rotations;
    end
    nElems   = analysis.getTotalElemsNumber();
    taskDim  = analysis.getTaskDim();
    nSupports   = nnz(analysis.supports);
    nLoadedDofs = nnz(analysis.Pnodal);

    if k == 1
        referenceElems     = model.mesh.elems;
        referenceModel     = model;
        referenceElemCount = nElems;
        referenceDofs      = analysis.ndofs;
        referenceTaskDim   = taskDim;
        fprintf('  Reference element count : %d\n', referenceElemCount);
        fprintf('  Reference task DOFs     : %d\n', referenceTaskDim);
    else
        assert(nElems == referenceElemCount, ...
            'Element-count mismatch in %s: got %d, expected %d.', ...
            cfg.name, nElems, referenceElemCount);
        assert(taskDim == referenceTaskDim, ...
            'DOF-count mismatch in %s: got %d, expected %d.', ...
            cfg.name, taskDim, referenceTaskDim);
        assertLinkedArmLayoutCompatible(model, referenceModel, cfg.name);
        assert(isequal(analysis.ndofs, referenceDofs), ...
            'DOF labels/order differ in configuration %s.', cfg.name);
    end

    models{k}   = model;
    analyses{k} = analysis;

    row.configName             = string(cfg.name);
    row.configLabel            = string(cfg.label);
    row.nNodes                 = size(model.mesh.nodes, 1);
    row.nElems                 = nElems;
    row.nTaskDofs              = taskDim;
    row.nSupportedDofs         = nSupports;
    row.nLoadedDofsBeforeSolve = nLoadedDofs;
    row.analysisClass          = string(class(analysis));
    row.useSecondOrderBuckling  = useSecondOrderBuckling;
    row.bucklingLambdaFactor    = bucklingLambdaFactor;
    row.sameConnectivityAsFirst = k == 1 || isequal(model.mesh.elems, referenceElems);
    row.sameLinkedLayoutAsFirst = true;
    row.sameDofsAsFirst         = true;
    setupRows{k, 1} = row;

    fprintf('  nodes=%d, elems=%d, taskDOFs=%d, supportedDOFs=%d, loadedDOFs=%d\n', ...
        row.nNodes, row.nElems, row.nTaskDofs, row.nSupportedDofs, row.nLoadedDofsBeforeSolve);
end

setupTable = struct2table(vertcat(setupRows{:}));
writetable(setupTable, fullfile(resultRoot, 'configuration_setup.csv'));

%% ---- Linked design space ---------------------------------------------------
H       = models{1}.halfSegmentNelems;
nElems  = referenceElemCount;
nCopies = nElems / H;

fprintf('\nLinked design space: H=%d reference elements, nElems=%d, nCopies=%d\n', ...
    H, nElems, nCopies);

assert(mod(nCopies, 2) == 0, ...
    'nElems/H = %d is not even; segmentToArm requires an even number of copies.', nCopies);

const_elems = armConstRingElementIds(models{1}, arm, "linked");
const_elems_full = armConstRingElementIds(models{1}, arm, "full");
fprintf('Const ring elems in linked rho: %d (full-arm copies: %d, end=%d, middle=%d)\n', ...
    numel(const_elems), numel(const_elems_full), arm.constEndRing, arm.constMiddleRing);

%% ---- Build and verify linked sensitivity map --------------------------------
fprintf('\nBuilding and verifying linked sensitivity map...\n');
map = buildLinkedSensitivityMap(H, nElems);

e_ref_test  = ceil(H / 2);
rho_ones    = ones(H, 1);
x_base      = models{1}.segmentToArm(rho_ones);
rhoPert     = rho_ones;
rhoPert(e_ref_test) = rhoPert(e_ref_test) + 1e-6;
x_pert      = models{1}.segmentToArm(rhoPert);
changedElems   = find(abs(x_pert - x_base) > 1e-10);
expectedAll    = sort([map.normalMap(:, e_ref_test); map.flippedMap(:, e_ref_test)]);
assert(isequal(changedElems, expectedAll), ...
    'Linked map verification failed for e_ref=%d.', e_ref_test);
fprintf('  Map verification PASSED: e_ref=%d changes exactly %d elements\n', ...
    e_ref_test, numel(expectedAll));

%% ---- Create stress-intensity topology optimization object -------------------
% The constructor builds the spatial filter matrix (O(nElems^2)) — may be slow
% for large meshes. The same FE analyses are passed; topology is shared.
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
topOpt.setConstElems(const_elems_full);
fprintf('  Done. Filter weights: [%d x %d] sparse.\n', ...
    size(topOpt.weights, 1), size(topOpt.weights, 2));
fprintf('  Aggregation class: %s\n', class(topOpt));

%% ---- Custom linked ESO loop ------------------------------------------------
fprintf('\nRunning linked ESO for up to %d iterations (target vf=%.3f)...\n', maxIter, VolFrac);

rho        = ones(H, 1);        % start from full material
erased_rho = false(H, 1);
rho(const_elems) = 1.0;

rhoHistory    = zeros(H, maxIter);
volHistory    = nan(maxIter, 1);
erasedHistory = nan(maxIter, 1);
maxStressHistory = nan(maxIter, nConfigs);
maxDisplacementHistory = nan(maxIter, nConfigs);
iterTimeHistory = nan(maxIter, 1);

nIter = 0;
for iter = 1:maxIter
    iterTic = tic;

    % -- 1. Expand rho -> x_arm -------------------------------------------------
    x_arm = models{1}.segmentToArm(rho);
    topOpt.x = x_arm;

    % -- 2. Solve all configs, compute filtered max-over-configs HM intensity ---
    nStressBefore = numel(topOpt.maxstress);
    nDisplacementBefore = numel(topOpt.maxdisplacement);
    ais_full = topOpt.weights * topOpt.computeAverageIntensities();
    % maxstress and maxdisplacement grew by nConfigs entries; grab latest values
    maxStressHistory(iter, :) = topOpt.maxstress(nStressBefore + 1 : end);
    maxDisplacementHistory(iter, :) = topOpt.maxdisplacement(nDisplacementBefore + 1 : end);

    % -- 3. Pull back to rho-space (mean over copies) --------------------------
    ais_rho = pullbackFullArmAverageIntensity(ais_full, H, nElems);

    % -- 4. Remove least-stressed rho elements (ESO rule) ----------------------
    notErasedID = setdiff(find(~erased_rho), const_elems);
    if ~isempty(notErasedID)
        ais_ne      = ais_rho(notErasedID);
        ais_range   = max(ais_ne) - min(ais_ne);
        removeList  = ais_ne < min(ais_ne) + ais_range * maxais;
        maxRemove   = round(H * max_elem_removal_factor);
        if sum(removeList) > maxRemove
            [~, ai]    = sort(ais_ne);
            removeList = false(numel(notErasedID), 1);
            removeList(ai(1:maxRemove)) = true;
        end
        erased_rho(notErasedID(removeList)) = true;
    end

    % -- 5. Update rho ----------------------------------------------------------
    rho(erased_rho)  = min(1, max(xminValue, ...
        rho(erased_rho) .* ais_rho(erased_rho) .^ penal));
    rho(~erased_rho) = 1;
    rho(const_elems) = 1.0;

    vol = mean(rho);
    iterTimeSec = toc(iterTic);

    % -- 6. Store history -------------------------------------------------------
    rhoHistory(:, iter)  = rho;
    volHistory(iter)     = vol;
    erasedHistory(iter)  = sum(erased_rho);
    iterTimeHistory(iter) = iterTimeSec;
    nIter                = iter;

    fprintf('%4d  vol=%.4f  erased=%d/%d  time=%.2fs  maxStress=[', ...
        iter, vol, sum(erased_rho), H, iterTimeSec);
    fprintf('%.3e ', maxStressHistory(iter, :));
    fprintf(']\n');

    % -- 7. Stop when target volume fraction reached ---------------------------
    if vol <= VolFrac
        fprintf('  Target volume fraction %.3f reached at iteration %d.\n', VolFrac, iter);
        break;
    end
end

%% ---- Final state and pass checks -------------------------------------------
rhoFinal            = rho;
x_arm_final         = models{1}.segmentToArm(rhoFinal);
finalVolumeFraction = mean(rhoFinal);
nErasedFinal        = sum(erased_rho);

volumeReached    = finalVolumeFraction <= VolFrac + 1e-2;
rhoNonuniform    = std(rhoFinal) > 0.03 && (max(rhoFinal) - min(rhoFinal)) > 0.15;
someErased       = nErasedFinal > 0;

fprintf('\nFinal checks\n');
fprintf('  Volume fraction reached   : %d (vf=%.6f, target=%.3f)\n', ...
    volumeReached, finalVolumeFraction, VolFrac);
fprintf('  Rho nonuniform            : %d (std=%.4f, range=%.4f)\n', ...
    rhoNonuniform, std(rhoFinal), max(rhoFinal) - min(rhoFinal));
fprintf('  Elements erased in rho    : %d / %d (%.1f%%)\n', ...
    nErasedFinal, H, 100 * nErasedFinal / H);

%% ---- Location-density diagnostics on full-arm topology ---------------------
locationStats = computeLocationDensityStats(models{1}, x_arm_final);
writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

%% ---- Save results ----------------------------------------------------------
% Trim history arrays to actual number of completed iterations
rhoHistory       = rhoHistory(:, 1:nIter);
volHistory       = volHistory(1:nIter);
erasedHistory    = erasedHistory(1:nIter);
maxStressHistory = maxStressHistory(1:nIter, :);
maxDisplacementHistory = maxDisplacementHistory(1:nIter, :);
iterTimeHistory  = iterTimeHistory(1:nIter);

history.iteration    = (1:nIter)';
history.volumeFraction = volHistory;
history.nErased      = erasedHistory;
history.iterationTimeSec = iterTimeHistory;
history.maxStress    = maxStressHistory;
history.maxDisplacement = maxDisplacementHistory;
history.x            = rhoHistory;
history.configNames  = string(cellfun(@(s) s.name,  configs, 'UniformOutput', false));
history.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

save(fullfile(resultRoot, 'result.mat'), ...
    'rhoFinal', 'x_arm_final', 'finalVolumeFraction', 'nErasedFinal', ...
    'erased_rho', 'history', 'configs', 'map', 'const_elems', 'const_elems_full', ...
    'H', 'nElems', 'nCopies', 'nIter', ...
    'VolFrac', 'penal', 'Rmin', 'maxais', 'maxIter', 'xminValue', ...
    'useParallel', 'parallelAvailable', 'useSecondOrderBuckling', 'bucklingLambdaFactor', ...
    'stressAggregation', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', 'Pz', ...
    'locationStats', 'rhoNonuniform', 'volumeReached', 'someErased');

% History CSV
histRows = struct();
for i = 1:nIter
    histRows(i).iteration     = i;
    histRows(i).volumeFraction = volHistory(i);
    histRows(i).nErased       = erasedHistory(i);
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
ppOpts.VolFrac      = VolFrac;
ppOpts.expand       = @(z) models{1}.segmentToArm(z);
ppOpts.penal        = penal;
ppOpts.useParallel  = useParallel;
ppOpts.resultRoot   = resultRoot;
ppOpts.configNames  = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
ppOpts.plotModels   = models;
ppOpts.plotConfigs  = configs;
ppOpts.runReanalysis = false;

postResult = postprocessStressIntensityResult(history, models{1}, analyses, ppOpts);

plotFinalTopology(models{1}, x_arm_final, resultRoot, postResult, analyses);
plotTopologyConfigurations(models, x_arm_final > 0.5, resultRoot, ...
    "final_threshold_rho_gt_05_by_config", "Final rho > 0.5 by configuration", configs);
plotLocationDensity(locationStats, resultRoot);

%% ---- Structural performance metrics ----------------------------------------
fprintf('\nComputing structural performance metrics...\n');
x_ref_C = ones(nElems, 1);
metrics_ref_C   = evaluateStructuralPerformance(analyses, x_ref_C,    1,     useParallel);
metrics_final_C = evaluateStructuralPerformance(analyses, x_arm_final, penal, useParallel);
saveStructuralMetricsCsv(metrics_ref_C, metrics_final_C, configs, finalVolumeFraction, resultRoot);

% Volume fraction convergence
fig = figure('Name', 'Test C volume history');
plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.2);
yline(VolFrac, '--r', sprintf('Target %.2f', VolFrac));
grid on;
xlabel('Iteration'); ylabel('Volume fraction');
title('Test C: linked ESO volume fraction');
saveas(fig, fullfile(resultRoot, 'volume_history.png'));
savefig(fig, fullfile(resultRoot, 'volume_history.fig'));
close(fig);

% Max stress per config
fig = figure('Name', 'Test C max stress history');
plot(history.iteration, history.maxStress, '-o', 'LineWidth', 1.0);
grid on;
xlabel('Iteration'); ylabel('Max HM stress');
legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), 'Location', 'best');
title('Test C: per-config max stress vs iteration');
saveas(fig, fullfile(resultRoot, 'stress_history.png'));
savefig(fig, fullfile(resultRoot, 'stress_history.fig'));
close(fig);

% Max displacement per config
fig = figure('Name', 'Test C max displacement history');
plot(history.iteration, history.maxDisplacement, '-o', 'LineWidth', 1.0);
grid on;
xlabel('Iteration'); ylabel('Max displacement magnitude');
legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), 'Location', 'best');
title('Test C: per-config max displacement vs iteration');
saveas(fig, fullfile(resultRoot, 'displacement_history.png'));
savefig(fig, fullfile(resultRoot, 'displacement_history.fig'));
close(fig);

% Reference rho bar plot
fig = figure('Name', 'Test C reference rho');
bar(rhoFinal, 'FaceColor', [0.72 0.33 0.28]);
xlabel('Reference half-segment element index');
ylabel('Density \rho');
title(sprintf('Test C: reference rho (H=%d, vf=%.3f)', H, mean(rhoFinal)));
ylim([0 1.05]);
exportgraphics(fig, fullfile(resultRoot, 'rho_final.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'rho_final.fig'));
close(fig);

%% ---- Optional: compare to testB (linked MMA) and testA (unlinked MMA) -----
testBPath = fullfile(scriptDir, 'results', 'testB_fullArmLinkedSIMP', 'result.mat');
testAPath = fullfile(scriptDir, 'results', 'testA_fullArmUnlinkedSIMP', 'result.mat');
if exist(testBPath, 'file') || exist(testAPath, 'file')
    fprintf('\nComparison\n');
    compRow.testCVolFrac  = finalVolumeFraction;
    compRow.testCnIter    = nIter;
    compRow.testCnErased  = nErasedFinal;
    if exist(testBPath, 'file')
        resB = load(testBPath, 'finalJ', 'finalVolumeFraction');
        fprintf('  Test B (linked MMA)   vf=%.6f  J=%.6e\n', ...
            resB.finalVolumeFraction, resB.finalJ);
        compRow.testBFinalJ   = resB.finalJ;
        compRow.testBVolFrac  = resB.finalVolumeFraction;
    end
    if exist(testAPath, 'file')
        resA = load(testAPath, 'finalJ', 'finalVolumeFraction');
        fprintf('  Test A (unlinked MMA) vf=%.6f  J=%.6e\n', ...
            resA.finalVolumeFraction, resA.finalJ);
        compRow.testAFinalJ  = resA.finalJ;
        compRow.testAVolFrac = resA.finalVolumeFraction;
    end
    writetable(struct2table(compRow), fullfile(resultRoot, 'comparison.csv'));
end

fprintf('\nSaved Test C outputs to %s\n', resultRoot);

%% ---- Pass criteria warnings ------------------------------------------------
if ~volumeReached
    warning('TestC:VolumeTargetNotReached', ...
        'Volume fraction %.6f did not reach target %.3f.', finalVolumeFraction, VolFrac);
end
if ~rhoNonuniform
    warning('TestC:RhoTooUniform', ...
        'Reference rho is too uniform after ESO: std=%.4f, range=%.4f.', ...
        std(rhoFinal), max(rhoFinal) - min(rhoFinal));
end
if ~someErased
    warning('TestC:NoElementsErased', ...
        'No rho elements were erased -- ESO loop did not run or remove anything.');
end

% testB_fullArmLinkedSIMP
% Test B: full-arm linked multi-configuration SIMP for modular Arm-Z design.
%
% Linking mode: rotation-aware (use_offset=false, nCircDiv=arm.nCircDiv).
%   resCirc is snapped to the nearest multiple of nCircDiv so every joint
%   angle is an exact integer multiple of 2*pi/resCirc.  Each half-segment
%   mesh is generated in its LOCAL frame (no phase trick); junction nodes
%   align automatically.  rho(e) then corresponds to the SAME local
%   circumferential position in every segment — a true repeating module.
%
% Formulation:
%   Design variables: rho [H x 1] — reference half-segment densities.
%   Full-arm expansion: x_arm = model.segmentToArm(rho)
%   Pattern: [rho; flip(rho); rho; flip(rho); ...]
%   (flip handles the reversed element ordering of 2b half-segments)
%
%   All configurations share the same x_arm; compliance and sensitivity are
%   computed on the full solid FEM for each configuration, then pulled back
%   to rho-space via the exact chain rule of the segmentToArm map.
%
% Objective:
%   J = (sum_k w_k * (C_k/C0_k)^pAgg)^(1/pAgg)
%
% Constraint:
%   mean(rho) <= VolFrac
%   (equivalent to mean(x_arm) <= VolFrac since mean is preserved by segmentToArm)
%
% Configurations (3):
%   1. Max M_z bending
%   2. Max M_s torsion
%   3. Max T_y shear

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(22, 'twister');

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

VolFrac = 0.40;
penal = 3.0;
pAgg = 4.0;
maxIter = 200;
xminValue = arm.mma.xminValue;
Rfilter = arm.Rfilter;
useParallel = license('test', 'Distrib_Computing_Toolbox');

fdElemCount = 5;
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

resultRoot = fullfile(scriptDir, 'results', 'testB_fullArmLinkedSIMP');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test B full-arm linked SIMP (modular Arm-Z)\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf('  VolFrac=%.3f, penal=%.2f, pAgg=%.2f, maxIter=%d, changeTol=%.1e, xmin=%.3f, parallel=%d\n', ...
    VolFrac, penal, pAgg, maxIter, changeTol, xminValue, useParallel);

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

    model    = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
        cfg.betas, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
    analysis = model.analysis;
    nElems   = analysis.getTotalElemsNumber();
    taskDim  = analysis.getTaskDim();
    nSupports    = nnz(analysis.supports);
    nLoadedDofs  = nnz(analysis.Pnodal);

    if k == 1
        referenceElems     = model.mesh.elems;
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
        assert(isequal(model.mesh.elems, referenceElems), ...
            'Mesh connectivity differs in configuration %s.', cfg.name);
        assert(isequal(analysis.ndofs, referenceDofs), ...
            'DOF labels/order differ in configuration %s.', cfg.name);
    end

    models{k}   = model;
    analyses{k} = analysis;

    row.configName          = string(cfg.name);
    row.configLabel         = string(cfg.label);
    row.nNodes              = size(model.mesh.nodes, 1);
    row.nElems              = nElems;
    row.nTaskDofs           = taskDim;
    row.nSupportedDofs      = nSupports;
    row.nLoadedDofsBeforeSolve = nLoadedDofs;
    row.sameConnectivityAsFirst = true;
    row.sameDofsAsFirst         = true;
    setupRows{k, 1}         = row;

    fprintf('  nodes=%d, elems=%d, taskDOFs=%d, supportedDOFs=%d, loadedDOFs=%d\n', ...
        row.nNodes, row.nElems, row.nTaskDofs, row.nSupportedDofs, row.nLoadedDofsBeforeSolve);
end

setupTable = struct2table(vertcat(setupRows{:}));
writetable(setupTable, fullfile(resultRoot, 'configuration_setup.csv'));

%% ---- Linked design space ---------------------------------------------------
H      = models{1}.halfSegmentNelems;
nElems = referenceElemCount;
nCopies = nElems / H;

fprintf('\nLinked design space: H=%d reference elements, nElems=%d, nCopies=%d\n', ...
    H, nElems, nCopies);

assert(mod(nCopies, 2) == 0, ...
    'nElems/H = %d is not even; segmentToArm requires an even number of copies.', nCopies);

xmin = xminValue * ones(H, 1);
xmax = ones(H, 1);
rho  = VolFrac * ones(H, 1);
const_elems = armConstRingElementIds(models{1}, arm, "linked");
xmin(const_elems) = 1.0;
xmax(const_elems) = 1.0;
rho(const_elems) = 1.0;
rho  = enforceVolumeFraction(rho, VolFrac, xmin, xmax);
fprintf('Const ring elems in linked rho: %d (end=%d, middle=%d)\n', ...
    numel(const_elems), arm.constEndRing, arm.constMiddleRing);

fprintf('\nBuilding linked reference sensitivity filter: Rfilter=%.6g, H=%d\n', ...
    Rfilter, H);
Wfilter = buildElementFilterMatrix(models{1}.mesh.nodes, models{1}.mesh.elems, ...
    (1:H)', Rfilter);

%% ---- Build and verify linked sensitivity map --------------------------------
fprintf('\nBuilding and verifying linked sensitivity map...\n');
map = buildLinkedSensitivityMap(H, nElems);

fprintf('  nArms=%d, nCopies=%d, normalMap size=[%s], flippedMap size=[%s]\n', ...
    map.nArms, map.nCopies, ...
    num2str(size(map.normalMap)), num2str(size(map.flippedMap)));

e_ref_test = ceil(H / 2);
x_base = models{1}.segmentToArm(rho);
rhoPert = rho;
rhoPert(e_ref_test) = rhoPert(e_ref_test) + 1e-6;
x_pert = models{1}.segmentToArm(rhoPert);
changedElems   = find(abs(x_pert - x_base) > 1e-10);
expectedNormal  = map.normalMap(:,  e_ref_test);
expectedFlipped = map.flippedMap(:, e_ref_test);
expectedAll     = sort([expectedNormal(:); expectedFlipped(:)]);
assert(isequal(changedElems, expectedAll), ...
    ['Linked map verification failed for e_ref=%d: ' ...
     '%d changed elements do not match %d expected.'], ...
    e_ref_test, numel(changedElems), numel(expectedAll));
fprintf('  Map verification PASSED: e_ref=%d changes exactly %d elements (%d normal + %d flipped)\n', ...
    e_ref_test, numel(expectedAll), numel(expectedNormal), numel(expectedFlipped));

%% ---- Compliance normalization and FD gradient check in rho-space -----------
fprintf('\nComputing initial compliance normalization at rho=VolFrac...\n');
x_arm = models{1}.segmentToArm(rho);
[J0, dJdx0, C0, Cinit] = evaluateObjectiveAndGradient(analyses, x_arm, penal, pAgg, weights, [], useParallel);
dJdrho0 = pullbackFullArmSensitivity(dJdx0, H, nElems);
fprintf('  Initial J = %.8e\n', J0);
for k = 1:nConfigs
    fprintf('    %-12s C0 = %.8e\n', configs{k}.name, C0(k));
end

fprintf('\nFinite-difference gradient test in rho-space (%d elements, h=%.1e)...\n', ...
    fdElemCount, fdStep);
fdLinked = finiteDifferenceGradientTestLinked(models{1}, analyses, rho, penal, pAgg, weights, ...
    C0, dJdrho0, fdElemCount, fdStep, xmin, xmax, useParallel);
fprintf('  FD (rho-space) max relative error = %.3e\n', fdLinked.maxRelativeError);
assert(fdLinked.maxRelativeError < fdRelTol, ...
    'Rho-space FD gradient check failed: max relative error %.3e exceeds %.3e.', ...
    fdLinked.maxRelativeError, fdRelTol);

fdTable = struct2table(fdLinked.rows);
writetable(fdTable, fullfile(resultRoot, 'finite_difference_gradient.csv'));

%% ---- MMA optimization on rho -----------------------------------------------
fprintf('\nRunning projected MMA on rho (%d variables) for %d iterations...\n', H, maxIter);

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
opts.expand = @(z) models{1}.segmentToArm(z);
opts.pullback = @(dx) pullbackFullArmSensitivity(dx, H, nElems);
opts.fixedDesignVariables = const_elems;
opts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
opts.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

optResult = solveSIMPComplianceVolumeMMA(analyses, rho, xmin, xmax, opts);

%% ---- Final state and pass checks -------------------------------------------
rhoFinal            = optResult.zFinal;
x_arm_final         = optResult.xFinal;
finalJ              = optResult.finalJ;
finalC              = optResult.finalC;
finalVolumeFraction = optResult.finalVolumeFraction;
finalChange         = optResult.finalChange;
nIter               = optResult.nIter;
history             = optResult.history;
objectiveWindowConverged = optResult.objectiveWindowConverged;
objectiveDecreased       = optResult.objectiveDecreased;
mmaAcceptable            = optResult.mmaAcceptable;
volumeActive             = abs(finalVolumeFraction - VolFrac) < 5.0e-3;

% Nonuniformity is assessed on rho (the actual design variable).
% By construction of segmentToArm, all half-segment location means equal
% mean(rho), so the testA location-diversity metric is zero for any linked
% design; instead we check that rho itself has structural variation.
rhoNonuniform        = std(rhoFinal) > 0.03 && (max(rhoFinal) - min(rhoFinal)) > 0.15;
topologyNonuniform   = rhoNonuniform;
locationDensityRange = max(rhoFinal) - min(rhoFinal);
locationDensityStd   = std(rhoFinal);

% Reuse topologyDiffersByArmLocation flag to signal internal rho variation
topologyDiffersByArmLocation = locationDensityRange > 0.03;

%% ---- Location-density diagnostics on full-arm topology ---------------------
locationStats = computeLocationDensityStats(models{1}, x_arm_final);
writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

fprintf('\nFinal checks\n');
fprintf('  Objective decreased          : %d (J0=%.6e, Jf=%.6e)\n', ...
    objectiveDecreased, history.J(1), finalJ);
fprintf('  Objective window converged   : %d\n', objectiveWindowConverged);
fprintf('  MMA acceptable               : %d\n', mmaAcceptable);
fprintf('  Volume constraint active     : %d (vf=%.6f)\n', volumeActive, finalVolumeFraction);
fprintf('  Rho nonuniform               : %d (std=%.4f, range=%.4f)\n', ...
    rhoNonuniform, std(rhoFinal), max(rhoFinal) - min(rhoFinal));

%% ---- Save histories and plots ----------------------------------------------

save(fullfile(resultRoot, 'result.mat'), ...
    'rhoFinal', 'x_arm_final', 'finalJ', 'finalC', ...
    'finalVolumeFraction', 'finalChange', ...
    'history', 'configs', 'map', ...
    'H', 'nElems', 'nCopies', 'nIter', ...
    'VolFrac', 'penal', 'pAgg', 'maxIter', 'minIter', 'changeTol', 'xminValue', 'Rfilter', 'useParallel', ...
    'const_elems', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', 'Pz', ...
    'locationStats', 'locationDensityRange', 'locationDensityStd', ...
    'topologyDiffersByArmLocation', 'objectiveDecreased', ...
    'objectiveWindowConverged', 'mmaAcceptable', 'volumeActive', ...
    'topologyNonuniform', 'rhoNonuniform');

saveHistoryCsv(resultRoot, history, configs);
saveSummaryCsv(resultRoot, fdLinked, history, rhoFinal, finalC, VolFrac, ...
    objectiveDecreased, objectiveWindowConverged, mmaAcceptable, volumeActive, ...
    topologyNonuniform, topologyDiffersByArmLocation, ...
    locationDensityRange, locationDensityStd);

%% ---- Post-processing: topology extraction ----------------------------------
fprintf('\nPost-processing topology extraction...\n');
ppOpts = struct();
ppOpts.penal       = penal;
ppOpts.pAgg        = pAgg;
ppOpts.weights     = weights;
ppOpts.VolFrac     = VolFrac;
ppOpts.fixedVars   = const_elems;
ppOpts.useParallel = useParallel;
ppOpts.resultRoot  = resultRoot;
ppOpts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
% Heaviside sharpening is only valid for projection-based SIMP (solveSIMPVolumeStressMMA).
% solveSIMPComplianceVolumeMMA uses a sensitivity filter only: z IS the physical density,
% so Wfilter*z has no density-sharpening meaning and collapses volume.
ppOpts.skipHeaviside = true;
% Set runReanalysis=true to enable expensive binary FE validation (~30 FEM solves).
ppOpts.runReanalysis = false;
% Set runSweep=true to generate a J-vs-V Pareto curve (~50 FEM solves).
ppOpts.runSweep = false;

postResult = postprocessSIMPResult(optResult, models{1}, analyses, Wfilter, ppOpts);

% Plot full-arm expanded topology with comparison panel
plotFinalTopology(models{1}, x_arm_final, resultRoot, postResult, analyses);
plotHistory(history, configs, resultRoot);
plotLocationDensity(locationStats, resultRoot);

%% ---- Structural performance metrics ----------------------------------------
fprintf('\nComputing structural performance metrics...\n');
x_ref_B = ones(nElems, 1);
metrics_ref_B   = evaluateStructuralPerformance(analyses, x_ref_B,    1,     useParallel);
metrics_final_B = evaluateStructuralPerformance(analyses, x_arm_final, penal, useParallel);
saveStructuralMetricsCsv(metrics_ref_B, metrics_final_B, configs, finalVolumeFraction, resultRoot);

% Plot reference module rho
figure('Visible', 'off');
bar(rhoFinal, 'FaceColor', [0.28 0.45 0.72]);
xlabel('Reference half-segment element index');
ylabel('Density \rho');
title(sprintf('Test B: reference half-segment density (H=%d, vf=%.3f)', ...
    H, mean(rhoFinal)));
ylim([0 1.05]);
exportgraphics(gcf, fullfile(resultRoot, 'rho_final.png'), 'Resolution', 200);
savefig(gcf, fullfile(resultRoot, 'rho_final.fig'));
close(gcf);

fprintf('\nSaved Test B outputs to %s\n', resultRoot);

%% ---- Optional: compare to Test A (unlinked) --------------------------------
testAPath = fullfile(scriptDir, 'results', 'testA_fullArmUnlinkedSIMP', 'result.mat');
if exist(testAPath, 'file')
    fprintf('\nComparison: Test B (linked) vs Test A (unlinked)\n');
    resA = load(testAPath, 'finalJ', 'finalVolumeFraction', 'VolFrac', 'penal', 'pAgg');
    fprintf('  Test A  final J = %.8e  vf = %.6f\n', resA.finalJ, resA.finalVolumeFraction);
    fprintf('  Test B  final J = %.8e  vf = %.6f\n', finalJ,      finalVolumeFraction);
    fprintf('  Ratio B/A       = %.4f  (>1 = cost of modularity constraint)\n', ...
        finalJ / resA.finalJ);
    compRow.testAFinalJ       = resA.finalJ;
    compRow.testAVolFrac      = resA.finalVolumeFraction;
    compRow.testBFinalJ       = finalJ;
    compRow.testBVolFrac      = finalVolumeFraction;
    compRow.ratioBoverA       = finalJ / resA.finalJ;
    compRow.linkedDesignVars  = H;
    compRow.unlinkedDesignVars= nElems;
    writetable(struct2table(compRow), fullfile(resultRoot, 'comparison_testA.csv'));
else
    fprintf('\n(Test A result not found; skipping A/B comparison.)\n');
end

%% ---- Pass criteria assertions ----------------------------------------------
assert(mmaAcceptable, ...
    'MMA did not converge by objective window and did not decrease J over %d iterations.', maxIter);
assert(volumeActive, ...
    'Volume constraint is not active enough: vf=%.6f, target=%.6f.', ...
    finalVolumeFraction, VolFrac);
assert(topologyNonuniform, ...
    'Reference module rho is too uniform: std=%.4f, range=%.4f.', ...
    std(rhoFinal), max(rhoFinal) - min(rhoFinal));

% testB_fullArmLinkedSIMP
% Test B: full-arm linked multi-configuration SIMP for modular Arm-Z design.
%
% Formulation:
%   Design variables: rho [H x 1] — reference half-segment densities.
%   Full-arm expansion: x_arm = model.segmentToArm(rho)
%   Pattern: [rho; flip(rho); rho; flip(rho); ...]
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

fdElemCount = 5;
fdStep = 1.0e-5;
fdRelTol = 1.0e-3;

moveLimit = 0.05;
minMoveLimit = 0.003;
moveDecay = 0.95;
mmaDamping = 0.25;
objectiveTol = 5.0e-3;

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

resultRoot = fullfile(scriptDir, 'results', 'testB_fullArmLinkedSIMP');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test B full-arm linked SIMP (modular Arm-Z)\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf('  VolFrac=%.3f, penal=%.2f, pAgg=%.2f, maxIter=%d, xmin=%.3f\n', ...
    VolFrac, penal, pAgg, maxIter, xminValue);

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
        cfg.betas, ShapeFn, true, Pz);
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
rho  = enforceVolumeFraction(rho, VolFrac, xmin, xmax);

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
[J0, dJdx0, C0, Cinit] = evaluateObjectiveAndGradient(analyses, x_arm, penal, pAgg, weights, []);
dJdrho0 = pullbackFullArmSensitivity(dJdx0, H, nElems);
fprintf('  Initial J = %.8e\n', J0);
for k = 1:nConfigs
    fprintf('    %-12s C0 = %.8e\n', configs{k}.name, C0(k));
end

fprintf('\nFinite-difference gradient test in rho-space (%d elements, h=%.1e)...\n', ...
    fdElemCount, fdStep);
fdLinked = finiteDifferenceGradientTestLinked(models{1}, analyses, rho, penal, pAgg, weights, ...
    C0, dJdrho0, fdElemCount, fdStep, xmin, xmax);
fprintf('  FD (rho-space) max relative error = %.3e\n', fdLinked.maxRelativeError);
assert(fdLinked.maxRelativeError < fdRelTol, ...
    'Rho-space FD gradient check failed: max relative error %.3e exceeds %.3e.', ...
    fdLinked.maxRelativeError, fdRelTol);

fdTable = struct2table(fdLinked.rows);
writetable(fdTable, fullfile(resultRoot, 'finite_difference_gradient.csv'));

%% ---- MMA optimization on rho -----------------------------------------------
fprintf('\nRunning projected MMA on rho (%d variables) for %d iterations...\n', H, maxIter);

rhoHistory = zeros(H, maxIter + 1);
rhoHistory(:, 1) = rho;
JHistory      = nan(maxIter + 1, 1);
JHistory(1)   = J0;
CHistory      = nan(maxIter + 1, nConfigs);
CHistory(1,:) = Cinit(:)';
volHistory    = nan(maxIter + 1, 1);
volHistory(1) = mean(rho);
changeHistory = nan(maxIter + 1, 1);
changeHistory(1) = 0;

m     = 1;
n     = H;
xold1 = rho;
xold2 = rho;
low   = zeros(n, 1);
upp   = ones(n, 1);
a0    = 1;
a     = 0;
c_mma = 1000;
d     = 0;
objectiveScale = 1.0 / max(abs(J0), eps);

for iter = 1:maxIter
    x_arm = models{1}.segmentToArm(rho);
    [J, dJdx, ~, ~] = evaluateObjectiveAndGradient(analyses, x_arm, penal, pAgg, weights, C0);
    dJdrho = pullbackFullArmSensitivity(dJdx, H, nElems);

    constr     = sum(rho) / (VolFrac * H) - 1.0;
    gradConstr = ones(1, H) / (VolFrac * H);

    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = mmasub2( ...
        m, n, iter, rho, xmin, xmax, xold1, xold2, ...
        objectiveScale * J,      objectiveScale * dJdrho,  0 * dJdrho, ...
        constr, gradConstr, 0 * gradConstr, ...
        low, upp, a0, a, c_mma, d);

    if iter > 1
        xold2 = xold1;
    end
    xold1 = rho;

    currentMoveLimit = max(minMoveLimit, moveLimit * moveDecay^(iter - 1));
    rhoCandidate = min(max(xmma, rho - currentMoveLimit), rho + currentMoveLimit);
    rhoCandidate = rho + mmaDamping * (rhoCandidate - rho);
    rhoCandidate = enforceVolumeFraction(rhoCandidate, VolFrac, xmin, xmax);

    change = max(abs(rhoCandidate - rho));
    rho    = rhoCandidate;

    x_arm = models{1}.segmentToArm(rho);
    [Jnew, Cnew] = evaluateObjectiveOnly(analyses, x_arm, penal, pAgg, weights, C0);

    rhoHistory(:, iter + 1) = rho;
    JHistory(iter + 1)      = Jnew;
    CHistory(iter + 1, :)   = Cnew(:)';
    volHistory(iter + 1)    = mean(rho);
    changeHistory(iter + 1) = change;

    fprintf('%4d  J=%12.6e  dJ/J0=% .3e  vf=%.4f  change=%.3e', ...
        iter, Jnew, (Jnew - JHistory(iter)) / max(abs(JHistory(iter)), eps), ...
        volHistory(iter + 1), change);
    for k = 1:nConfigs
        fprintf('  C_%s=%.3e', configs{k}.name, Cnew(k));
    end
    fprintf('\n');
end

%% ---- Final state and pass checks -------------------------------------------
rhoFinal            = rho;
x_arm_final         = models{1}.segmentToArm(rhoFinal);
finalJ              = JHistory(maxIter + 1);
finalC              = CHistory(maxIter + 1, :)';
finalVolumeFraction = mean(rhoFinal);
finalChange         = changeHistory(maxIter + 1);

objectiveWindowConverged = objectiveHistoryConverged(JHistory, objectiveTol);
objectiveDecreased       = finalJ < JHistory(1);
mmaAcceptable            = objectiveWindowConverged || objectiveDecreased;
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
    objectiveDecreased, JHistory(1), finalJ);
fprintf('  Objective window converged   : %d\n', objectiveWindowConverged);
fprintf('  MMA acceptable               : %d\n', mmaAcceptable);
fprintf('  Volume constraint active     : %d (vf=%.6f)\n', volumeActive, finalVolumeFraction);
fprintf('  Rho nonuniform               : %d (std=%.4f, range=%.4f)\n', ...
    rhoNonuniform, std(rhoFinal), max(rhoFinal) - min(rhoFinal));

%% ---- Save histories and plots ----------------------------------------------
history.iteration   = (0:maxIter)';
history.J           = JHistory;
history.C           = CHistory;
history.volumeFraction = volHistory;
history.change      = changeHistory;
history.x           = rhoHistory;   % design-variable history is rho (H x nIter+1)
history.configNames = string(cellfun(@(s) s.name,  configs, 'UniformOutput', false));
history.configLabels= string(cellfun(@(s) s.label, configs, 'UniformOutput', false));
history.weights     = weights;
history.C0          = C0;

save(fullfile(resultRoot, 'result.mat'), ...
    'rhoFinal', 'x_arm_final', 'finalJ', 'finalC', ...
    'finalVolumeFraction', 'finalChange', ...
    'history', 'configs', 'map', ...
    'H', 'nElems', 'nCopies', ...
    'VolFrac', 'penal', 'pAgg', 'maxIter', 'xminValue', ...
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

% Plot full-arm expanded topology (x_arm_final)
plotFinalTopology(models{1}, x_arm_final, resultRoot);
plotHistory(history, configs, resultRoot);
plotLocationDensity(locationStats, resultRoot);

% Plot reference module rho
figure('Visible', 'off');
bar(rhoFinal, 'FaceColor', [0.28 0.45 0.72]);
xlabel('Reference half-segment element index');
ylabel('Density \rho');
title(sprintf('Test B: reference half-segment density (H=%d, vf=%.3f)', ...
    H, mean(rhoFinal)));
ylim([0 1.05]);
exportgraphics(gcf, fullfile(resultRoot, 'rho_final.png'), 'Resolution', 200);
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

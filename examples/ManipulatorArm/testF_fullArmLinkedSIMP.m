% testF_fullArmLinkedSIMP
% Test F: full-arm linked inverse SIMP formulation for modular Arm-Z design.
%
% Objective:
%   minimize mean(rho)
%
% Constraints:
%   C_k <= C_coeff * C0_k for every load configuration k
%
% C0_k is the full-density linked-design compliance. At the minimum-volume
% solution, active compliance constraints behave as C_k = C_coeff*C0_k.

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(77, 'twister');

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
useParallel = license('test', 'Distrib_Computing_Toolbox');

C_coeff = 2.0;          % allow up to 2x full-density compliance
constraintTol = 5.0e-3;
changeTol = arm.mma.changeTol;
minIter = 10;

moveLimit = arm.mma.moveLimit;
minMoveLimit = arm.mma.minMoveLimit;
moveDecay = arm.mma.moveDecay;
mmaDamping = arm.mma.mmaDamping;

%configs = armLoadConfigs("sixPlusTension");
configs = armLoadConfigs("six");
nConfigs = numel(configs);

resultRoot = fullfile(scriptDir, 'results', 'testF_fullArmLinkedSIMP');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test F full-arm linked inverse SIMP\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf('  Objective=min volume, constraints C_k <= %.3f*C0_k, maxIter=%d, parallel=%d\n', ...
    C_coeff, maxIter, useParallel);

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

%% ---- Initial compliance targets --------------------------------------------
fprintf('\nComputing full-density linked compliance targets...\n');
x_arm = models{1}.segmentToArm(rho);
[C0, ~] = evaluateComplianceSet(analyses, x_arm, penal, useParallel);
Ctarget = C_coeff * max(C0, eps);
for k = 1:nConfigs
    fprintf('  %-12s C0=%.8e  Ctarget=%.8e\n', configs{k}.name, C0(k), Ctarget(k));
end

%% ---- MMA optimization on rho ------------------------------------------------
fprintf('\nRunning inverse projected MMA on rho (%d variables) for up to %d iterations...\n', ...
    H, maxIter);

opts = struct();
opts.penal = penal;
opts.C_coeff = C_coeff;
opts.maxIter = maxIter;
opts.minIter = minIter;
opts.changeTol = changeTol;
opts.constraintTol = constraintTol;
opts.moveLimit = moveLimit;
opts.minMoveLimit = minMoveLimit;
opts.moveDecay = moveDecay;
opts.mmaDamping = mmaDamping;
opts.useParallel = useParallel;
opts.sensitivityFilter = Wfilter;
opts.expand = @(z) models{1}.segmentToArm(z);
opts.pullback = @(dx) pullbackComplianceSet(dx, H, nElems);
opts.fixedDesignVariables = const_elems;
opts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
opts.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

optResult = solveSIMPVolumeComplianceMMA(analyses, rho, xmin, xmax, opts);

%% ---- Final state and checks -------------------------------------------------
rhoFinal = optResult.zFinal;
x_arm_final = optResult.xFinal;
finalVolumeFraction = optResult.finalVolumeFraction;
finalC = optResult.finalC;
finalConstraint = optResult.finalConstraint;
finalChange = optResult.finalChange;
constraintsSatisfied = optResult.constraintsSatisfied;
nIter = optResult.nIter;
history = optResult.history;
C0 = optResult.C0;
Ctarget = optResult.Ctarget;
rhoNonuniform = std(rhoFinal) > 0.03 && (max(rhoFinal) - min(rhoFinal)) > 0.15;
topologyNonuniform = rhoNonuniform;
locationDensityRange = max(rhoFinal) - min(rhoFinal);
locationDensityStd = std(rhoFinal);
topologyDiffersByArmLocation = locationDensityRange > 0.03;

locationStats = computeLocationDensityStats(models{1}, x_arm_final);
writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

fprintf('\nFinal checks\n');
fprintf('  Final volume fraction        : %.6f\n', finalVolumeFraction);
fprintf('  Compliance constraints ok    : %d (max g=%.3e)\n', constraintsSatisfied, max(finalConstraint));
fprintf('  Rho nonuniform               : %d (std=%.4f, range=%.4f)\n', ...
    rhoNonuniform, std(rhoFinal), max(rhoFinal) - min(rhoFinal));

%% ---- Save histories and plots ----------------------------------------------
save(fullfile(resultRoot, 'result.mat'), ...
    'rhoFinal', 'x_arm_final', 'finalC', 'finalVolumeFraction', 'finalConstraint', 'finalChange', ...
    'history', 'configs', 'map', 'H', 'nElems', 'nCopies', 'nIter', ...
    'C0', 'Ctarget', 'C_coeff', ...
    'penal', 'maxIter', 'minIter', 'changeTol', 'constraintTol', 'xminValue', 'Rfilter', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', 'useParallel', 'Pz', ...
    'const_elems', ...
    'locationStats', 'locationDensityRange', 'locationDensityStd', ...
    'topologyDiffersByArmLocation', 'constraintsSatisfied', 'topologyNonuniform', 'rhoNonuniform');

saveInverseHistoryCsv(resultRoot, history, configs);
saveInverseSummaryCsv(resultRoot, history, rhoFinal, finalC, finalConstraint, C_coeff, ...
    constraintsSatisfied, topologyNonuniform, topologyDiffersByArmLocation, ...
    locationDensityRange, locationDensityStd);

plotTopology(models{1}, x_arm_final, resultRoot, struct('smoothed', true, 'saveFig', true));
exportTopology(models{1}, x_arm_final, resultRoot);
plotTopologyConfigurations(models, x_arm_final > 0.5, resultRoot, ...
    "final_threshold_rho_gt_05_by_config", "Final rho > 0.5 by configuration", configs);
plotInverseHistory(history, configs, resultRoot, 'Test F');
plotLocationDensity(locationStats, resultRoot);

fig = figure('Name', 'Test F reference rho');
bar(rhoFinal, 'FaceColor', [0.28 0.45 0.72]);
xlabel('Reference half-segment element index');
ylabel('Density \rho');
title(sprintf('Test F: reference half-segment density (H=%d, vf=%.3f)', ...
    H, mean(rhoFinal)));
ylim([0 1.05]);
exportgraphics(fig, fullfile(resultRoot, 'rho_final.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'rho_final.fig'));
close(fig);

fprintf('\nSaved Test F outputs to %s\n', resultRoot);

if ~constraintsSatisfied
    warning('TestF:ComplianceConstraintsViolated', ...
        'Final compliance constraints are violated: max g=%.3e.', max(finalConstraint));
end

%% ---- Local helpers ----------------------------------------------------------
function dCdrho = pullbackComplianceSet(dCdx, H, nElems)
    nConfigs = size(dCdx, 2);
    dCdrho = zeros(H, nConfigs);
    for k = 1:nConfigs
        dCdrho(:, k) = pullbackFullArmSensitivity(dCdx(:, k), H, nElems);
    end
end

function [C, dC] = evaluateComplianceSet(analyses, x, penal, useParallel)
    nConfigs = numel(analyses);
    C = zeros(nConfigs, 1);
    dC = zeros(numel(x), nConfigs);
    useParallel = useParallel && license('test', 'Distrib_Computing_Toolbox');
    if useParallel
        parfor k = 1:nConfigs
            [C_k, dC_k] = computeComplianceAndGradient(analyses{k}, x, penal);
            C(k) = C_k;
            dC(:, k) = dC_k;
        end
    else
        for k = 1:nConfigs
            [C(k), dC(:, k)] = computeComplianceAndGradient(analyses{k}, x, penal);
        end
    end
end

function saveInverseHistoryCsv(resultRoot, history, configs)
    nRows = numel(history.iteration);
    rows = repmat(struct(), nRows, 1);
    for i = 1:nRows
        rows(i).iteration = history.iteration(i);
        rows(i).volumeFraction = history.volumeFraction(i);
        rows(i).change = history.change(i);
        rows(i).iterationTimeSec = history.iterationTimeSec(i);
        for k = 1:numel(configs)
            rows(i).(['C_' configs{k}.name]) = history.C(i, k);
            rows(i).(['g_' configs{k}.name]) = history.constraint(i, k);
        end
    end
    writetable(struct2table(rows), fullfile(resultRoot, 'history.csv'));
end

function saveInverseSummaryCsv(resultRoot, history, finalX, finalC, finalConstraint, C_coeff, ...
        constraintsSatisfied, topologyNonuniform, topologyDiffersByArmLocation, ...
        locationDensityRange, locationDensityStd)
    row.status = "ok";
    row.nElements = numel(finalX);
    row.nIterations = numel(history.iteration) - 1;
    row.C_coeff = C_coeff;
    row.finalVolumeFraction = mean(finalX);
    row.finalChange = history.change(end);
    row.maxConstraint = max(finalConstraint);
    row.constraintsSatisfied = constraintsSatisfied;
    row.topologyNonuniform = topologyNonuniform;
    row.finalDensityStd = std(finalX);
    row.finalDensityRange = max(finalX) - min(finalX);
    row.rhoGt05 = mean(finalX > 0.5);
    row.rhoGt07 = mean(finalX > 0.7);
    row.finalCMax = max(finalC);
    row.topologyDiffersByArmLocation = topologyDiffersByArmLocation;
    row.locationMeanDensityRange = locationDensityRange;
    row.locationMeanDensityStd = locationDensityStd;
    writetable(struct2table(row), fullfile(resultRoot, 'summary.csv'));
end

function plotInverseHistory(history, configs, resultRoot, testName)
    fig = figure('Name', [testName ' volume history']);
    plot(history.iteration, history.volumeFraction, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('Volume fraction');
    title([testName ': volume objective']);
    saveas(fig, fullfile(resultRoot, 'volume_history.png'));
    savefig(fig, fullfile(resultRoot, 'volume_history.fig'));
    close(fig);

    fig = figure('Name', [testName ' compliance history']);
    plot(history.iteration, history.C, '-o', 'LineWidth', 1.0);
    grid on; xlabel('Iteration'); ylabel('Compliance');
    legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), 'Location', 'best');
    title([testName ': compliance constraints']);
    saveas(fig, fullfile(resultRoot, 'compliance_history.png'));
    savefig(fig, fullfile(resultRoot, 'compliance_history.fig'));
    close(fig);

    fig = figure('Name', [testName ' constraint history']);
    plot(history.iteration, history.constraint, '-o', 'LineWidth', 1.0);
    yline(0, '--k', 'Constraint');
    grid on; xlabel('Iteration'); ylabel('C/Ctarget - 1');
    legend(cellfun(@(s) s.label, configs, 'UniformOutput', false), 'Location', 'best');
    title([testName ': normalized compliance constraints']);
    saveas(fig, fullfile(resultRoot, 'constraint_history.png'));
    savefig(fig, fullfile(resultRoot, 'constraint_history.fig'));
    close(fig);
end

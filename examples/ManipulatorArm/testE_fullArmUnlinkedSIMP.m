% testE_fullArmUnlinkedSIMP
% Test E: full-arm unlinked inverse SIMP formulation.
%
% Objective:
%   minimize mean(x)
%
% Constraints:
%   C_k <= C_coeff * C0_k for every load configuration k
%
% C0_k is the full-density initial compliance. At the minimum-volume
% solution, active compliance constraints behave as C_k = C_coeff*C0_k.

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(55, 'twister');

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

resultRoot = fullfile(scriptDir, 'results', 'testE_fullArmUnlinkedSIMP');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

fprintf('Test E full-arm unlinked inverse SIMP\n');
fprintf('  Result root: %s\n', resultRoot);
fprintf('  Objective=min volume, constraints C_k <= %.3f*C0_k, maxIter=%d, parallel=%d\n', ...
    C_coeff, maxIter, useParallel);

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
        cfg.betas, ShapeFn, true, Pz, arm.constEndRing, arm.constMiddleRing);
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

%% ---- Initial compliance targets --------------------------------------------
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

fprintf('\nComputing full-density compliance targets...\n');
[C0, ~] = evaluateComplianceSet(analyses, x, penal, useParallel);
Ctarget = C_coeff * max(C0, eps);
for k = 1:nConfigs
    fprintf('  %-12s C0=%.8e  Ctarget=%.8e\n', configs{k}.name, C0(k), Ctarget(k));
end

%% ---- MMA optimization -------------------------------------------------------
fprintf('\nRunning inverse projected MMA on x (%d variables) for up to %d iterations...\n', ...
    nDesign, maxIter);

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
opts.fixedDesignVariables = const_elems;
opts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
opts.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));

optResult = solveSIMPVolumeComplianceMMA(analyses, x, xmin, xmax, opts);

%% ---- Final state and checks -------------------------------------------------
finalX = optResult.zFinal;
finalVolumeFraction = optResult.finalVolumeFraction;
finalC = optResult.finalC;
finalConstraint = optResult.finalConstraint;
finalChange = optResult.finalChange;
constraintsSatisfied = optResult.constraintsSatisfied;
nIter = optResult.nIter;
history = optResult.history;
C0 = optResult.C0;
Ctarget = optResult.Ctarget;
topologyNonuniform = std(finalX) > 0.03 && (max(finalX) - min(finalX)) > 0.15;

locationStats = computeLocationDensityStats(models{1}, finalX);
locationDensityRange = max([locationStats.meanDensity]) - min([locationStats.meanDensity]);
locationDensityStd = std([locationStats.meanDensity]);
topologyDiffersByArmLocation = locationDensityRange > 0.03;
writetable(struct2table(locationStats), fullfile(resultRoot, 'location_density.csv'));

fprintf('\nFinal checks\n');
fprintf('  Final volume fraction        : %.6f\n', finalVolumeFraction);
fprintf('  Compliance constraints ok    : %d (max g=%.3e)\n', constraintsSatisfied, max(finalConstraint));
fprintf('  Topology nonuniform          : %d (std=%.4f, range=%.4f)\n', ...
    topologyNonuniform, std(finalX), max(finalX) - min(finalX));

%% ---- Save histories and plots ----------------------------------------------
save(fullfile(resultRoot, 'result.mat'), ...
    'finalX', 'finalC', 'finalVolumeFraction', 'finalConstraint', 'finalChange', ...
    'history', 'configs', 'nIter', 'C0', 'Ctarget', 'C_coeff', ...
    'penal', 'maxIter', 'minIter', 'changeTol', 'constraintTol', 'xminValue', 'Rfilter', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', 'useParallel', 'Pz', ...
    'const_elems', ...
    'locationStats', 'locationDensityRange', 'locationDensityStd', ...
    'topologyDiffersByArmLocation', 'constraintsSatisfied', 'topologyNonuniform');

saveInverseHistoryCsv(resultRoot, history, configs);
saveInverseSummaryCsv(resultRoot, history, finalX, finalC, finalConstraint, C_coeff, ...
    constraintsSatisfied, topologyNonuniform, topologyDiffersByArmLocation, ...
    locationDensityRange, locationDensityStd);

plotFinalTopology(models{1}, finalX, resultRoot, [], analyses);
plotInverseHistory(history, configs, resultRoot, 'Test E');
plotLocationDensity(locationStats, resultRoot);

fprintf('\nSaved Test E outputs to %s\n', resultRoot);

if ~constraintsSatisfied
    warning('TestE:ComplianceConstraintsViolated', ...
        'Final compliance constraints are violated: max g=%.3e.', max(finalConstraint));
end

%% ---- Local helpers ----------------------------------------------------------
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

function testCurveParamSegmentOptimization(varargin)
% testCurveParamSegmentOptimization
% Curve-parametrized linked-segment density optimization for the Arm-Z model.
%
% The density field is generated on one reference half-segment from helical,
% axial, and circumferential curve families, then expanded to the full arm
% through ManipulatorModel3D.segmentToArm.  This is a low-dimensional,
% manufacturability-oriented alternative to element-wise SIMP.
%
% Objective modes:
%   "stress"            minimize stress aggregate at fixed VolFrac  (fminsearch)
%   "uz"                minimize tip displacement at fixed VolFrac   (fminsearch)
%   "combined"          weighted stress + displacement               (fminsearch)
%   "min_volume_stress" minimize VF s.t. stress  <= stressConstraintRatio * full_pipe  (fmincon)
%   "min_volume_disp"   minimize VF s.t. tipUz   <= dispConstraintRatio  * full_pipe  (fmincon)
%   "min_volume"        minimize VF s.t. both stress and disp constraints              (fmincon)
%
% Usage examples:
%   testCurveParamSegmentOptimization('objectiveMode','min_volume_stress','stressConstraintRatio',3.0)
%   testCurveParamSegmentOptimization('objectiveMode','min_volume','stressConstraintRatio',3.0,'dispConstraintRatio',8.0)

close all; clc;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(31, 'twister');
runOpts = parseRunOptions(varargin{:});

%% ---- Geometry and optimizer settings ---------------------------------------
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

objectiveMode = runOpts.objectiveMode;
isMinVolumeMode = startsWith(objectiveMode, "min_volume");
stressWeight = 0.70;
dispWeight = 0.30;
pElem = 12.0;
pConfigStress = 4.0;
pConfigDisp = 4.0;

maxIter = runOpts.maxIter;
maxFunEvals = runOpts.maxFunEvals;

configs = armLoadConfigs("statistical");
if isfinite(runOpts.configLimit)
    configs = configs(1:min(numel(configs), runOpts.configLimit));
end
nConfigs = numel(configs);

if isempty(gcp('nocreate'))
    parpool('local', min(nConfigs, feature('numcores')));
end

% Auto-generate result tag for min_volume modes if not manually set.
if strlength(runOpts.resultTag) == 0 && isMinVolumeMode
    tag = strrep(objectiveMode, 'min_volume_', '') ;
    if tag == "min_volume", tag = "vol"; end
    if isfinite(runOpts.stressConstraintRatio)
        tag = tag + sprintf("_s%.1f", runOpts.stressConstraintRatio);
    end
    if isfinite(runOpts.dispConstraintRatio)
        tag = tag + sprintf("_d%.1f", runOpts.dispConstraintRatio);
    end
    runOpts.resultTag = tag;
end

resultName = 'testCurveParamSegmentOptimization';
if strlength(runOpts.resultTag) > 0
    resultName = resultName + "_" + runOpts.resultTag;
end
resultRoot = fullfile(scriptDir, 'results', char(resultName));
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

fprintf('Curve-parametrized linked segment optimization\n');
fprintf('  objective=%s, penal=%.2f, pElem=%.1f\n', objectiveMode, penal, pElem);
if isMinVolumeMode
    fprintf('  stressConstraintRatio=%.2f, dispConstraintRatio=%.2f\n', ...
        runOpts.stressConstraintRatio, runOpts.dispConstraintRatio);
else
    fprintf('  VolFrac=%.3f\n', VolFrac);
end
fprintf('  Result root: %s\n', resultRoot);

%% ---- Build full-arm configurations -----------------------------------------
models = cell(nConfigs, 1);
analyses = cell(nConfigs, 1);
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
        referenceElemCount = nElems;
        referenceDofs = analysis.ndofs;
        referenceTaskDim = taskDim;
        fprintf('  Reference element count: %d\n', referenceElemCount);
        fprintf('  Reference task DOFs    : %d\n', referenceTaskDim);
    else
        assert(nElems == referenceElemCount, 'Element-count mismatch in %s.', cfg.name);
        assert(taskDim == referenceTaskDim, 'DOF-count mismatch in %s.', cfg.name);
        assert(isequal(analysis.ndofs, referenceDofs), ...
            'DOF labels/order differ in configuration %s.', cfg.name);
        assertLinkedArmLayoutCompatible(model, models{1}, cfg.name);
    end

    models{k} = model;
    analyses{k} = analysis;
end

modelRef = models{1};
H = modelRef.halfSegmentNelems;
constRef = armConstRingElementIds(modelRef, arm, "linked");
constFull = armConstRingElementIds(modelRef, arm, "full");

fprintf('\nReference half-segment elements: %d\n', H);
fprintf('Constant linked ring elements : %d\n', numel(constRef));
fprintf('Constant full-arm ring elems  : %d\n', numel(constFull));

%% ---- Curve-density setup ---------------------------------------------------
curveOpts = struct();
curveOpts.VolFrac = VolFrac;
curveOpts.penal = penal;
curveOpts.constRefElems = constRef;
curveOpts.constFullElems = constFull;
curveOpts.rhoMin = arm.mma.xminValue;
curveOpts.elemSize = arm.nominalElementSize;
curveOpts.pElem = pElem;
curveOpts.pConfigStress = pConfigStress;
curveOpts.pConfigDisp = pConfigDisp;
curveOpts.objectiveMode = selectInnerObjectiveMode(objectiveMode);
curveOpts.stressWeight = stressWeight;
curveOpts.dispWeight = dispWeight;
curveOpts.resultRoot = resultRoot;
curveOpts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
curveOpts.configLabels = string(cellfun(@(s) s.label, configs, 'UniformOutput', false));
curveOpts.configWeights = buildConfigWeights(configs, runOpts);

paramBounds = defaultCurveParamBounds(arm);
p0 = defaultCurveParamInitial(arm);
u0 = curveParamsToUnconstrained(p0, paramBounds);

fprintf('\nInitial curve parameters:\n');
disp(struct2table(curveParamsToDisplayStruct(p0)));

%% ---- Reference metrics and optimization ------------------------------------
evalCount = 0;
historyRows = struct([]);

if isMinVolumeMode
    %% -- Option B: fmincon, minimize VF subject to stress/disp constraints --

    fprintf('\nComputing full-pipe reference metrics...\n');
    rhoFullPipe = ones(referenceElemCount, 1);
    [fullPipeMetrics, ~] = evaluateLinkedDensityMetrics(analyses, rhoFullPipe, curveOpts);
    stressFullPipe = max(fullPipeMetrics.maxHM);
    dispFullPipe   = max(abs(fullPipeMetrics.tipUz));

    stressLimit = runOpts.stressConstraintRatio * stressFullPipe;
    dispLimit   = runOpts.dispConstraintRatio   * dispFullPipe;
    fprintf('Full-pipe: maxHM=%.4e Pa, max|tipUz|=%.4e m\n', stressFullPipe, dispFullPipe);
    fprintf('Stress limit: %.4e Pa (x%.2f)   Disp limit: %.4e m (x%.2f)\n', ...
        stressLimit, runOpts.stressConstraintRatio, dispLimit, runOpts.dispConstraintRatio);

    % Cache to avoid double FEM evaluation when fmincon calls obj + nonlcon
    % at the same point (which it does every iteration).
    cachedTheta   = [];
    cachedMetrics = [];

    vfMin = runOpts.vfMin;
    % Auto-init: start at VF=1.0 (full pipe, feasible) so SQP reduces VF
    % from a feasible starting point rather than climbing out of infeasibility.
    if isempty(runOpts.vfInitial)
        vfInitial = 1.0;
    else
        vfInitial = runOpts.vfInitial;
    end
    fprintf('Starting VF: %.4f\n', vfInitial);
    theta0 = [u0(:); vfInitial];
    lb = [-Inf(numel(u0), 1); vfMin];
    ub = [ Inf(numel(u0), 1); 1.0 ];

    fminconOpts = optimoptions('fmincon', ...
        'Algorithm',                  'sqp', ...
        'Display',                    'iter', ...
        'MaxIterations',              maxIter, ...
        'MaxFunctionEvaluations',     maxFunEvals, ...
        'OptimalityTolerance',        1e-3, ...
        'ConstraintTolerance',        1e-3, ...
        'FiniteDifferenceStepSize',   0.05, ...
        'FiniteDifferenceType',       'central');

    fprintf('\nRunning fmincon (SQP) to minimize volume fraction...\n');
    if maxIter <= 0 || maxFunEvals <= 1
        fprintf('Skipping fmincon because MaxIter=%d and MaxFunEvals=%d.\n', maxIter, maxFunEvals);
        thetaBest = theta0;
    else
        [thetaBest, ~] = fmincon(@objectiveMinVol, theta0, [], [], [], [], lb, ub, ...
            @constraintMinVol, fminconOpts);
    end

    uBest   = thetaBest(1:end-1);
    vfBest  = thetaBest(end);
    pBest   = unconstrainedToCurveParams(uBest, paramBounds);
    JBest   = vfBest;

    curveOpts.VolFrac = vfBest;
    [rhoRefBest, rhoFullBest, bestDensityInfo] = buildCurveLinkedDensity(pBest, modelRef, curveOpts);
    [bestMetrics, bestRaw] = evaluateLinkedDensityMetrics(analyses, rhoFullBest, curveOpts);

    % Use full-pipe metrics as the comparison baseline for summary CSV.
    baselineMetrics = fullPipeMetrics;
    baselineRaw     = [];
    rhoUniformRef   = rhoFullPipe(1:H);
    rhoUniformFull  = rhoFullPipe;

else
    %% -- Existing modes: fminsearch at fixed VolFrac -------------------------

    fprintf('\nEvaluating uniform baseline at VolFrac=%.3f...\n', VolFrac);
    rhoUniformFull = uniformDensityWithFixedVolume(referenceElemCount, VolFrac, ...
        arm.mma.xminValue, constFull);
    rhoUniformRef = rhoUniformFull(1:H);
    rhoUniformFull(constFull) = 1.0;
    [baselineMetrics, baselineRaw] = evaluateLinkedDensityMetrics(analyses, rhoUniformFull, curveOpts);
    printCurveMetrics('Uniform baseline', baselineMetrics, configs);
    curveOpts.stressRef = baselineMetrics.stressAggregateByConfig;
    curveOpts.dispRef   = abs(baselineMetrics.tipUz);
    [baselineMetrics, baselineRaw] = evaluateLinkedDensityMetrics(analyses, rhoUniformFull, curveOpts);

    options = optimset('Display', 'iter', ...
        'MaxIter',     maxIter, ...
        'MaxFunEvals', maxFunEvals, ...
        'TolX',        1.0e-3, ...
        'TolFun',      1.0e-3);

    fprintf('\nRunning fminsearch over bounded curve parameters...\n');
    if maxIter <= 0 || maxFunEvals <= 1
        fprintf('Skipping fminsearch because MaxIter=%d and MaxFunEvals=%d.\n', maxIter, maxFunEvals);
        uBest = u0;
        JBest = objectiveFromUnconstrained(u0);
    else
        [uBest, JBest] = fminsearch(@objectiveFromUnconstrained, u0, options);
    end
    pBest = unconstrainedToCurveParams(uBest, paramBounds);
    [rhoRefBest, rhoFullBest, bestDensityInfo] = buildCurveLinkedDensity(pBest, modelRef, curveOpts);
    [bestMetrics, bestRaw] = evaluateLinkedDensityMetrics(analyses, rhoFullBest, curveOpts);

    vfBest = bestDensityInfo.fullVolumeFraction;
end

fprintf('\nBest curve parameters:\n');
disp(struct2table(curveParamsToDisplayStruct(pBest)));
fprintf('Best density volume: linked=%.5f, full=%.5f\n', ...
    bestDensityInfo.linkedVolumeFraction, bestDensityInfo.fullVolumeFraction);
printCurveMetrics('Best curve design', bestMetrics, configs);
if isMinVolumeMode
    achievedStressRatio = max(bestMetrics.maxHM) / stressFullPipe;
    achievedDispRatio   = max(abs(bestMetrics.tipUz)) / dispFullPipe;
    fprintf('Achieved: VF=%.4f, stress ratio=%.3f (limit %.2f), disp ratio=%.3f (limit %.2f)\n', ...
        vfBest, achievedStressRatio, runOpts.stressConstraintRatio, ...
        achievedDispRatio, runOpts.dispConstraintRatio);
end

%% ---- Save results ----------------------------------------------------------
saveVars = {'pBest', 'JBest', 'p0', 'vfBest', ...
    'rhoRefBest', 'rhoFullBest', 'bestMetrics', 'bestRaw', ...
    'baselineMetrics', 'baselineRaw', 'rhoUniformRef', 'rhoUniformFull', ...
    'historyRows', 'configs', 'curveOpts', 'paramBounds', 'arm'};
if isMinVolumeMode
    stressConstraintRatio = runOpts.stressConstraintRatio; %#ok<NASGU>
    dispConstraintRatio   = runOpts.dispConstraintRatio;   %#ok<NASGU>
    saveVars = [saveVars, {'stressConstraintRatio', 'dispConstraintRatio', ...
        'stressLimit', 'dispLimit', 'stressFullPipe', 'dispFullPipe', 'fullPipeMetrics'}];
end
save(fullfile(resultRoot, 'result.mat'), saveVars{:});

if ~isempty(historyRows)
    writetable(struct2table(historyRows), fullfile(resultRoot, 'history.csv'));
end
writeCurveSummaryCsv(resultRoot, pBest, bestMetrics, baselineMetrics, configs);

fig = figure('Color', 'white', 'Name', 'curve density reference half-segment');
plotElementDensityField(modelRef, [rhoRefBest; 0*rhoFullBest(H+1:end)]);
title('Curve-parametrized density on linked reference elements', 'Interpreter', 'none');
exportgraphics(fig, fullfile(resultRoot, 'rho_reference_curve.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'rho_reference_curve.fig'));

fig = figure('Color', 'white', 'Name', 'curve density full arm');
plotElementDensityField(modelRef, rhoFullBest);
title('Curve-parametrized linked density on full arm', 'Interpreter', 'none');
exportgraphics(fig, fullfile(resultRoot, 'rho_full_curve.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'rho_full_curve.fig'));

fprintf('\nSaved results to: %s\n', resultRoot);

postprocessCurveParamResult(resultRoot);

%% ---- Nested functions ------------------------------------------------------

function J = objectiveMinVol(theta)
    J = theta(end);  % minimize VF directly
end

function [c, ceq] = constraintMinVol(theta)
    % Return one inequality per configuration so SQP receives individual
    % gradients for every load case, not just the single worst-case max.
    m = evalForMinVol(theta);
    c = m.maxHM(:) / stressLimit - 1;          % [nConfigs x 1]
    if isfinite(dispLimit)
        c = [c; abs(m.tipUz(:)) / dispLimit - 1];  % append [nConfigs x 1]
    end
    ceq = [];
end

function m = evalForMinVol(theta)
    % Cache: rebuild only when theta changes (fmincon calls obj + nonlcon
    % separately but at the same point each iteration).
    if ~isequal(theta, cachedTheta)
        params = unconstrainedToCurveParams(theta(1:end-1), paramBounds);
        vf     = max(vfMin, min(1.0, theta(end)));
        localOpts          = curveOpts;
        localOpts.VolFrac  = vf;
        [~, rhoFull, densityInfo] = buildCurveLinkedDensity(params, modelRef, localOpts);
        [metrics, ~]       = evaluateLinkedDensityMetrics(analyses, rhoFull, localOpts);
        cachedTheta   = theta;
        cachedMetrics = metrics;

        evalCount = evalCount + 1;
        sRatio = max(metrics.maxHM)      / stressFullPipe;
        dRatio = max(abs(metrics.tipUz)) / dispFullPipe;
        fprintf('eval %4d: VF=%.4f, stress_ratio=%.3f (lim=%.2f), disp_ratio=%.3f (lim=%.2f)\n', ...
            evalCount, densityInfo.fullVolumeFraction, sRatio, runOpts.stressConstraintRatio, ...
            dRatio, runOpts.dispConstraintRatio);

        row = params;
        row.eval         = evalCount;
        row.VF           = densityInfo.fullVolumeFraction;
        row.stressRatio  = sRatio;
        row.dispRatio    = dRatio;
        row.maxHM        = max(metrics.maxHM);
        row.maxAbsUz     = max(abs(metrics.tipUz));
        historyRows = [historyRows; row]; %#ok<AGROW>
    end
    m = cachedMetrics;
end

function J = objectiveFromUnconstrained(u)
    evalCount = evalCount + 1;
    params = unconstrainedToCurveParams(u, paramBounds);
    [rhoRef, rhoFull, densityInfo] = buildCurveLinkedDensity(params, modelRef, curveOpts); %#ok<ASGLU>
    [metrics, raw] = evaluateLinkedDensityMetrics(analyses, rhoFull, curveOpts); %#ok<NASGU>
    J = metrics.objective;

    row = params;
    row.eval = evalCount;
    row.J = J;
    row.stressObjective = metrics.stressObjective;
    row.dispObjective = metrics.dispObjective;
    row.trueMaxHM = max(metrics.maxHM);
    row.maxAbsUz = max(abs(metrics.tipUz));
    row.fullVolumeFraction = densityInfo.fullVolumeFraction;
    for kk = 1:numel(configs)
        hmField = char("maxHM_" + string(configs{kk}.name));
        uzField = char("tipUz_" + string(configs{kk}.name));
        row.(hmField) = metrics.maxHM(kk);
        row.(uzField) = metrics.tipUz(kk);
    end
    historyRows = [historyRows; row]; %#ok<AGROW>

    fprintf(['eval %4d: J=%.6e, S=%.6e, U=%.6e, maxHM=%.6e, ' ...
        'max|uz|=%.6e, vf=%.4f\n'], ...
        evalCount, J, metrics.stressObjective, metrics.dispObjective, ...
        max(metrics.maxHM), max(abs(metrics.tipUz)), densityInfo.fullVolumeFraction);
end

function rho = uniformDensityWithFixedVolume(nElems, targetVf, rhoMin, fixedIds)
    rho = zeros(nElems, 1);
    fixed = false(nElems, 1);
    fixed(fixedIds) = true;
    free = ~fixed;
    freeValue = (targetVf * nElems - nnz(fixed)) / nnz(free);
    freeValue = min(1, max(rhoMin, freeValue));
    rho(free) = freeValue;
    rho(fixed) = 1.0;
end

function opts = parseRunOptions(varargin)
    opts = struct();
    opts.objectiveMode         = "stress";
    opts.maxIter               = 80;
    opts.maxFunEvals           = 400;
    opts.configLimit           = Inf;
    opts.torsionWeight         = 1.0;
    opts.bendingWeight         = 1.0;
    opts.shearWeight           = 1.0;
    opts.resultTag             = "";
    opts.stressConstraintRatio = Inf;
    opts.dispConstraintRatio   = Inf;
    opts.vfInitial             = [];   % auto: 1.0 for min_volume, unused otherwise
    opts.vfMin                 = 0.05;

    if numel(varargin) == 1 && strcmpi(string(varargin{1}), "quick")
        opts.maxIter = 0;
        opts.maxFunEvals = 1;
        opts.configLimit = 2;
        return;
    end

    assert(mod(numel(varargin), 2) == 0, ...
        'Use name-value options, or the single option "quick".');
    for ii = 1:2:numel(varargin)
        name = lower(char(string(varargin{ii})));
        value = varargin{ii + 1};
        switch name
            case "objectivemode"
                opts.objectiveMode = string(value);
            case "maxiter"
                opts.maxIter = value;
            case "maxfunevals"
                opts.maxFunEvals = value;
            case "configlimit"
                opts.configLimit = value;
            case "torsionweight"
                opts.torsionWeight = value;
            case "bendingweight"
                opts.bendingWeight = value;
            case "shearweight"
                opts.shearWeight = value;
            case "resulttag"
                opts.resultTag = string(value);
            case "stressconstraintratio"
                opts.stressConstraintRatio = value;
            case "dispconstraintratio"
                opts.dispConstraintRatio = value;
            case "vfinitial"
                opts.vfInitial = value;
            case "vfmin"
                opts.vfMin = value;
            otherwise
                error('Unknown option "%s".', name);
        end
    end

    % Validate constraints are set for min_volume modes.
    if startsWith(string(opts.objectiveMode), "min_volume")
        if isinf(opts.stressConstraintRatio) && isinf(opts.dispConstraintRatio)
            error(['min_volume modes require at least one constraint. ' ...
                'Set stressConstraintRatio and/or dispConstraintRatio.']);
        end
    end
end

function mode = selectInnerObjectiveMode(outerMode)
    % min_volume modes bypass metrics.objective; use "stress" so
    % evaluateLinkedDensityMetrics stays valid and stressRef/dispRef
    % fields are populated (even though they won't be used as objective).
    known = ["stress", "uz", "combined"];
    if any(outerMode == known)
        mode = outerMode;
    else
        mode = "stress";
    end
end

function weights = buildConfigWeights(configs, runOpts)
    weights = ones(numel(configs), 1);
    for jj = 1:numel(configs)
        name = string(configs{jj}.name);
        if contains(name, "torsion")
            weights(jj) = runOpts.torsionWeight;
        elseif contains(name, "bending")
            weights(jj) = runOpts.bendingWeight;
        elseif contains(name, "shear")
            weights(jj) = runOpts.shearWeight;
        end
    end
    if sum(weights) <= 0
        weights(:) = 1;
    end
    weights = weights / sum(weights);
end

function s = curveParamsToDisplayStruct(params)
    % Return a struct containing only the 14 optimisation parameters in
    % their canonical order, suppressing legacy fields like angleDeg.
    names = curveParamNames();
    s = struct();
    for i = 1:numel(names)
        s.(names{i}) = params.(names{i});
    end
end

end

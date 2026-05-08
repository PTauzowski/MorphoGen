function testCurveParamRobustBetaOptimization(varargin)
% testCurveParamRobustBetaOptimization
% Robust topology optimization: minimize volume fraction subject to stress
% and displacement constraints holding for all adversarially discovered
% joint-angle configurations (constraint generation over discrete beta grid).
%
% Usage examples:
%   testCurveParamRobustBetaOptimization()
%   testCurveParamRobustBetaOptimization('deltaDeg', 45, 'stressConstraintRatio', 3.0)
%   testCurveParamRobustBetaOptimization('deltaDeg', 90, 'maxCGIter', 5, 'topK', 5)

close all; clc;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(31, 'twister');
runOpts = parseRobustRunOptions(varargin{:});

%% ---- Geometry defaults --------------------------------------------------------
arm    = armModelDefaults("thin");
E      = arm.E;
nu     = arm.nu;
R      = arm.R;
r      = arm.r;
h_seg  = arm.h_seg;
alpha  = arm.alpha;
res    = arm.res;
res_th = arm.res_th;
Pz     = arm.Pz;
ShapeFn = arm.ShapeFn;

penal         = 3.0;
pElem         = 12.0;
pConfigStress = 4.0;
pConfigDisp   = 4.0;

%% ---- Result folder ------------------------------------------------------------
tag = sprintf('robust_d%g_s%.1f', runOpts.deltaDeg, runOpts.stressConstraintRatio);
if isfinite(runOpts.dispConstraintRatio)
    tag = [tag, sprintf('_d%.1f', runOpts.dispConstraintRatio)];
end
if strlength(runOpts.resultTag) > 0
    tag = char(runOpts.resultTag);
end
resultName = 'testCurveParamRobustBetaOptimization_' + string(tag);
resultRoot = fullfile(scriptDir, 'results', char(resultName));
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

fprintf('Robust beta optimization\n');
fprintf('  deltaDeg=%g  stressConstraintRatio=%.2f  dispConstraintRatio=%.2f\n', ...
    runOpts.deltaDeg, runOpts.stressConstraintRatio, runOpts.dispConstraintRatio);
fprintf('  maxCGIter=%d  topK=%d  propLevel=%d\n', ...
    runOpts.maxCGIter, runOpts.topK, runOpts.propLevel);
fprintf('  Result root: %s\n\n', resultRoot);

%% ---- Build nominal reference model and full-pipe metrics ----------------------
nominalBeta = zeros(1, 7);
modelRef = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
    nominalBeta, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
analysisRef = modelRef.analysis;

H        = modelRef.halfSegmentNelems;
constRef  = armConstRingElementIds(modelRef, arm, "linked");
constFull = armConstRingElementIds(modelRef, arm, "full");

fprintf('Reference half-segment elements: %d\n', H);
fprintf('Constant linked ring elements  : %d\n', numel(constRef));

%% ---- Curve-density options ----------------------------------------------------
curveOpts = struct();
curveOpts.penal          = penal;
curveOpts.constRefElems  = constRef;
curveOpts.constFullElems = constFull;
curveOpts.rhoMin         = arm.mma.xminValue;
curveOpts.elemSize       = arm.nominalElementSize;
curveOpts.pElem          = pElem;
curveOpts.pConfigStress  = pConfigStress;
curveOpts.pConfigDisp    = pConfigDisp;
curveOpts.objectiveMode  = 'stress';
curveOpts.resultRoot     = resultRoot;

paramBounds = defaultCurveParamBounds(arm);
p0          = defaultCurveParamInitial(arm);

%% ---- Full-pipe reference for constraint limits --------------------------------
if isempty(gcp('nocreate'))
    parpool('local', feature('numcores'));
end

fprintf('Computing full-pipe reference metrics...\n');
rhoFullPipe = ones(analysisRef.getTotalElemsNumber(), 1);
% Use reference analysis to get single-config full-pipe values.
[fpMetrics, ~] = evaluateLinkedDensityMetrics({analysisRef}, rhoFullPipe, curveOpts);
stressFullPipe = fpMetrics.maxHM(1);
dispFullPipe   = abs(fpMetrics.tipUz(1));

stressLimit = runOpts.stressConstraintRatio * stressFullPipe;
dispLimit   = runOpts.dispConstraintRatio   * dispFullPipe;

fprintf('Full-pipe: maxHM=%.4e Pa, |tipUz|=%.4e m\n', stressFullPipe, dispFullPipe);
fprintf('Stress limit: %.4e Pa  (x%.2f)\n', stressLimit, runOpts.stressConstraintRatio);
fprintf('Disp   limit: %.4e m   (x%.2f)\n\n', dispLimit, runOpts.dispConstraintRatio);

%% ---- Constraint-generation robust optimization --------------------------------
cgOpts.maxCGIter       = runOpts.maxCGIter;
cgOpts.topK            = runOpts.topK;
cgOpts.propLevel       = runOpts.propLevel;
cgOpts.penal           = penal;
cgOpts.verbose         = true;
cgOpts.innerSolverOpts.maxIter     = runOpts.maxInnerIter;
cgOpts.innerSolverOpts.maxFunEvals = runOpts.maxInnerFunEvals;
cgOpts.innerSolverOpts.vfInitial   = 1.0;
cgOpts.innerSolverOpts.vfMin       = 0.05;

[pBest, vfBest, adversarialConfigs, historyRows] = constraintGenerationRobust( ...
    modelRef, arm, curveOpts, stressLimit, dispLimit, ...
    paramBounds, p0, runOpts.deltaDeg, cgOpts);

%% ---- Evaluate final design on full active config set --------------------------
fprintf('Evaluating final design...\n');
curveOpts.VolFrac = vfBest;
[rhoRefBest, rhoFullBest, bestDensityInfo] = buildCurveLinkedDensity(pBest, modelRef, curveOpts);

% Build nominal analysis for final metric report.
[bestMetrics, bestRaw] = evaluateLinkedDensityMetrics({analysisRef}, rhoFullBest, curveOpts);

fprintf('\nBest curve parameters:\n');
disp(struct2table(curveParamsToDisplayStruct(pBest)));
fprintf('VF = %.4f (full arm: %.5f)\n', vfBest, bestDensityInfo.fullVolumeFraction);
fprintf('Nominal: s_ratio=%.3f  d_ratio=%.3f  (limit 1.00; >1 means violated)\n', ...
    bestMetrics.maxHM(1) / stressLimit, ...
    abs(bestMetrics.tipUz(1)) / max(dispLimit, eps));

%% ---- Save results -------------------------------------------------------------
saveVars = {'pBest', 'vfBest', 'p0', ...
    'rhoRefBest', 'rhoFullBest', 'bestMetrics', 'bestRaw', ...
    'adversarialConfigs', 'historyRows', 'curveOpts', 'paramBounds', 'arm', ...
    'stressLimit', 'dispLimit', 'stressFullPipe', 'dispFullPipe', ...
    'stressConstraintRatio', 'dispConstraintRatio'};
stressConstraintRatio = runOpts.stressConstraintRatio; %#ok<NASGU>
dispConstraintRatio   = runOpts.dispConstraintRatio;   %#ok<NASGU>
save(fullfile(resultRoot, 'result.mat'), saveVars{:});

if ~isempty(historyRows)
    writetable(struct2table(historyRows), fullfile(resultRoot, 'history.csv'));
end
if ~isempty(adversarialConfigs)
    saveAdversarialConfigSummary(resultRoot, adversarialConfigs, stressLimit, dispLimit);
end

%% ---- Plots --------------------------------------------------------------------
fig = figure('Color', 'white', 'Name', 'curve density reference half-segment');
plotElementDensityField(modelRef, [rhoRefBest; zeros(numel(rhoFullBest)-H, 1)]);
title('Robust curve-parametrized density (reference segment)', 'Interpreter', 'none');
exportgraphics(fig, fullfile(resultRoot, 'rho_reference_curve.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'rho_reference_curve.fig'));

fig = figure('Color', 'white', 'Name', 'curve density full arm');
plotElementDensityField(modelRef, rhoFullBest);
title('Robust curve-parametrized linked density (full arm)', 'Interpreter', 'none');
exportgraphics(fig, fullfile(resultRoot, 'rho_full_curve.png'), 'Resolution', 200);
savefig(fig, fullfile(resultRoot, 'rho_full_curve.fig'));

fprintf('\nSaved results to: %s\n', resultRoot);
end

% =========================================================================
% Local helpers
% =========================================================================

function opts = parseRobustRunOptions(varargin)
    opts.deltaDeg             = 90;
    opts.stressConstraintRatio = 4.0;
    opts.dispConstraintRatio   = Inf;
    opts.maxCGIter            = 10;
    opts.topK                 = 3;
    opts.propLevel            = 1;
    opts.maxInnerIter         = 80;
    opts.maxInnerFunEvals     = 400;
    opts.resultTag            = "";

    for k = 1:2:numel(varargin)
        key = char(varargin{k});
        val = varargin{k+1};
        switch key
            case 'deltaDeg',              opts.deltaDeg             = val;
            case 'stressConstraintRatio', opts.stressConstraintRatio = val;
            case 'dispConstraintRatio',   opts.dispConstraintRatio   = val;
            case 'maxCGIter',             opts.maxCGIter            = val;
            case 'topK',                  opts.topK                 = val;
            case 'propLevel',             opts.propLevel            = val;
            case 'maxInnerIter',          opts.maxInnerIter         = val;
            case 'maxInnerFunEvals',      opts.maxInnerFunEvals     = val;
            case 'resultTag',             opts.resultTag            = string(val);
            otherwise
                warning('testCurveParamRobustBetaOptimization: unknown option "%s"', key);
        end
    end

    if isinf(opts.dispConstraintRatio)
        opts.dispConstraintRatio = Inf;
    end
end

function saveAdversarialConfigSummary(resultRoot, adversarialConfigs, stressLimit, dispLimit)
% Write adversarial_config_summary.csv with all verified candidates.
    rows = struct([]);
    for ci = 1:numel(adversarialConfigs)
        c = adversarialConfigs(ci);
        row.ci          = ci;
        row.stressRatio = c.stressRatio;
        row.dispRatio   = c.dispRatio;
        row.violated    = double(c.violated);
        row.maxHM_MPa   = c.maxHM / 1e6;
        row.tipUz_m     = c.tipUz;
        betaVec = c.beta;
        for ji = 1:numel(betaVec)
            row.(sprintf('beta%d', ji)) = betaVec(ji);
        end
        rows = [rows; row]; %#ok<AGROW>
    end
    if ~isempty(rows)
        writetable(struct2table(rows), ...
            fullfile(resultRoot, 'adversarial_config_summary.csv'));
    end
end

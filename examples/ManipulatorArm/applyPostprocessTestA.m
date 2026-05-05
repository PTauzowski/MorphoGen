% applyPostprocessTestA
% Apply topology post-processing to an already-saved Test A result.
%
% Loads results/testA_fullArmUnlinkedSIMP/result.mat, rebuilds the model
% and filter matrix (fast), recomputes C0 (one FEM solve), then calls
% postprocessSIMPResult.  The full optimisation does NOT need to be re-run.

clear; close all; clc;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

resultRoot = fullfile(scriptDir, 'results', 'testA_fullArmUnlinkedSIMP');
matFile    = fullfile(resultRoot, 'result.mat');
assert(exist(matFile, 'file') == 2, 'result.mat not found at %s', matFile);

fprintf('Loading result from %s\n', matFile);
r = load(matFile, ...
    'finalX', 'finalJ', 'finalC', 'VolFrac', 'penal', 'pAgg', ...
    'Rfilter', 'const_elems', 'configs', ...
    'E', 'nu', 'R', 'r', 'h_seg', 'alpha', 'res', 'res_th', 'Pz');

%% ---- Rebuild model and analyses from saved geometry -----------------------
arm     = armModelDefaults("thin");
configs = r.configs;
nConfigs = numel(configs);
weights  = ones(nConfigs, 1) / nConfigs;

fprintf('Rebuilding %d configurations...\n', nConfigs);
models   = cell(nConfigs, 1);
analyses = cell(nConfigs, 1);

for k = 1:nConfigs
    cfg      = configs{k};
    model    = ManipulatorModel3D(r.E, r.nu, r.h_seg, r.R, r.r, r.res, r.res_th, ...
        r.alpha, cfg.betas, arm.ShapeFn, true, r.Pz, arm.constEndRing, arm.constMiddleRing);
    models{k}   = model;
    analyses{k} = model.analysis;
end

nElems = analyses{1}.getTotalElemsNumber();
fprintf('  %d elements per configuration.\n', nElems);

%% ---- Rebuild filter matrix ------------------------------------------------
fprintf('Building element filter matrix (Rfilter=%.4g, nElems=%d)...\n', r.Rfilter, nElems);
Wfilter = buildElementFilterMatrix(models{1}.mesh.nodes, models{1}.mesh.elems, ...
    (1:nElems)', r.Rfilter);
fprintf('  Done.\n');

%% ---- Recompute C0 at initial density -------------------------------------
xminValue = arm.mma.xminValue;
xmin = xminValue * ones(nElems, 1);
xmax = ones(nElems, 1);
xmin(r.const_elems) = 1.0;
xmax(r.const_elems) = 1.0;

x_init = r.VolFrac * ones(nElems, 1);
x_init(r.const_elems) = 1.0;
x_init = enforceVolumeFraction(x_init, r.VolFrac, xmin, xmax);

fprintf('Recomputing C0 at initial density (one FEM solve per config)...\n');
useParallel = license('test', 'Distrib_Computing_Toolbox');
[~, ~, C0, ~] = evaluateObjectiveAndGradient( ...
    analyses, x_init, r.penal, r.pAgg, weights, [], useParallel);
C0 = max(C0, eps);
fprintf('  C0 = [%s]\n', num2str(C0', '%.4e '));

%% ---- Assemble optResult-compatible struct ---------------------------------
optResult.zFinal = r.finalX;
optResult.xFinal = r.finalX;   % unlinked: z = x (no segmentToArm)
optResult.C0     = C0;

%% ---- Post-processing options ----------------------------------------------
ppOpts = struct();
ppOpts.penal        = r.penal;
ppOpts.pAgg         = r.pAgg;
ppOpts.weights      = weights;
ppOpts.VolFrac      = r.VolFrac;
ppOpts.fixedVars    = r.const_elems;
ppOpts.useParallel  = useParallel;
ppOpts.resultRoot   = resultRoot;
ppOpts.configNames  = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
ppOpts.skipHeaviside = true;  % sensitivity-filtered SIMP: z IS physical density
ppOpts.runReanalysis = false;
ppOpts.runSweep     = false;

%% ---- Run post-processing --------------------------------------------------
fprintf('\nRunning topology post-processing...\n');
postResult = postprocessSIMPResult(optResult, models{1}, analyses, Wfilter, ppOpts);

plotFinalTopology(models{1}, r.finalX, resultRoot, postResult, analyses);

fprintf('\nPost-processing complete.  Outputs saved to:\n  %s\n', resultRoot);
fprintf('Best method : %s  (J_simp=%.6f, V=%.4f)\n', ...
    postResult.best.label, postResult.best.J_simp, postResult.best.volFrac);

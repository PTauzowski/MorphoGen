% repaintTestBResults
% Rebuild Test B models from result.mat and regenerate plot artifacts only.

clear; close all; clc;
clear classes;

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

resultRoot = fullfile(scriptDir, 'results', 'testB_fullArmLinkedSIMP');
resultPath = fullfile(resultRoot, 'result.mat');
assert(exist(resultPath, 'file') == 2, 'Missing Test B result.mat: %s', resultPath);

fprintf('Repainting Test B plots from %s\n', resultPath);
r = load(resultPath);

% Remove stale plot files only. Keep result.mat and CSV outputs.
plotPatterns = ["final_*.png", "final_*.fig", ...
                "postprocess_*.png", "postprocess_*.fig", ...
                "convergence_history.png", "convergence_history.fig", ...
                "location_density.png", "location_density.fig", ...
                "rho_final.png", "rho_final.fig"];
for p = plotPatterns
    files = dir(fullfile(resultRoot, p));
    for i = 1:numel(files)
        delete(fullfile(files(i).folder, files(i).name));
    end
end

arm = armModelDefaults("thin");
configs = r.configs;
nConfigs = numel(configs);

models = cell(nConfigs, 1);
analyses = cell(nConfigs, 1);
referenceModel = [];
referenceDofs = [];

fprintf('Rebuilding %d linked configurations for plotting...\n', nConfigs);
for k = 1:nConfigs
    cfg = configs{k};
    fprintf('  %d/%d %s betas=%s\n', k, nConfigs, cfg.name, mat2str(cfg.betas));
    model = ManipulatorModel3D(r.E, r.nu, r.h_seg, r.R, r.r, r.res, r.res_th, r.alpha, ...
        cfg.betas, arm.ShapeFn, false, r.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
    analysis = model.analysis;

    if k == 1
        referenceModel = model;
        referenceDofs = analysis.ndofs;
    else
        assertLinkedArmLayoutCompatible(model, referenceModel, cfg.name);
        assert(isequal(analysis.ndofs, referenceDofs), ...
            'DOF labels/order differ in configuration %s.', cfg.name);
    end

    models{k} = model;
    analyses{k} = analysis;
end

H = models{1}.halfSegmentNelems;
nElems = analyses{1}.getTotalElemsNumber();
assert(numel(r.rhoFinal) == H, 'rhoFinal length does not match H.');
assert(numel(r.x_arm_final) == nElems, 'x_arm_final length does not match model element count.');

fprintf('Rebuilding linked sensitivity filter...\n');
Wfilter = buildElementFilterMatrix(models{1}.mesh.nodes, models{1}.mesh.elems, ...
    (1:H)', r.Rfilter);

fprintf('Recomputing postprocess candidates and plots...\n');
optResult = struct();
optResult.zFinal = r.rhoFinal;
optResult.xFinal = r.x_arm_final;
if isfield(r.history, 'C') && size(r.history.C, 1) >= 1
    optResult.C0 = r.history.C(1, :)';
else
    optResult.C0 = max(r.finalC(:), eps);
end

ppOpts = struct();
ppOpts.penal       = r.penal;
ppOpts.pAgg        = r.pAgg;
ppOpts.weights     = r.history.weights(:);
ppOpts.VolFrac     = r.VolFrac;
ppOpts.fixedVars   = r.const_elems;
ppOpts.useParallel = false;
ppOpts.resultRoot  = resultRoot;
ppOpts.configNames = string(cellfun(@(s) s.name, configs, 'UniformOutput', false));
ppOpts.plotModels  = models;
ppOpts.plotConfigs = configs;
ppOpts.skipHeaviside = true;
ppOpts.runReanalysis = false;
ppOpts.runSweep = false;

postResult = postprocessSIMPResult(optResult, models{1}, analyses, Wfilter, ppOpts);

plotFinalTopology(models{1}, r.x_arm_final, resultRoot, postResult, []);
plotTopologyConfigurations(models, r.x_arm_final > 0.5, resultRoot, ...
    "final_threshold_rho_gt_05_by_config", "Final rho > 0.5 by configuration", configs);
plotHistory(r.history, configs, resultRoot);
plotLocationDensity(r.locationStats, resultRoot);

fig = figure('Visible', 'off', 'Name', 'Test B reference rho');
bar(r.rhoFinal, 'FaceColor', [0.28 0.45 0.72]);
xlabel('Reference half-segment element index');
ylabel('Density \rho');
title(sprintf('Test B: reference half-segment density (H=%d, vf=%.3f)', ...
    H, mean(r.rhoFinal)));
ylim([0 1.05]);
exportgraphics(fig, fullfile(resultRoot, 'rho_final.png'), 'Resolution', 200);
close(fig);

fprintf('Repainted Test B plots in %s\n', resultRoot);

function postprocessCurveParamResult(resultRoot)
% postprocessCurveParamResult  Reanalyse and export curve-parametric designs.
%
% Saves 3D topology files (STL / coloured OBJ) and metric CSVs.
% Curve-family fields are exported for one module only (single segment).

    close all; clc;

    scriptDir = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(scriptDir, '..', '..');
    addpath(genpath(projectRoot));

    if nargin < 1 || isempty(resultRoot)
        resultRoot = fullfile(scriptDir, 'results', 'testCurveParamSegmentOptimization');
    end
    resultRoot = char(string(resultRoot));
    resultFile = fullfile(resultRoot, 'result.mat');
    assert(exist(resultFile, 'file') == 2, 'Missing result file: %s', resultFile);

    S = load(resultFile);
    isRobustResult = isfield(S, 'adversarialConfigs') && ~isempty(S.adversarialConfigs);
    if isRobustResult
        required = {'pBest', 'rhoFullBest', 'curveOpts', 'arm'};
    else
        required = {'pBest', 'rhoFullBest', 'configs', 'curveOpts', 'arm'};
    end
    for i = 1:numel(required)
        assert(isfield(S, required{i}), 'result.mat lacks "%s".', required{i});
    end

    pBest       = S.pBest;
    rhoFullBest = S.rhoFullBest(:);
    curveOpts   = S.curveOpts;
    arm         = S.arm;
    if isRobustResult
        configs = {struct('name', 'nominal', 'label', 'Nominal (straight)', ...
                          'betas', zeros(1, 7))};
    else
        configs = S.configs;
    end
    adversarialConfigs = [];
    if isfield(S, 'adversarialConfigs')
        adversarialConfigs = S.adversarialConfigs;
    end
    stressLimit = localOpt(S, 'stressLimit', Inf);
    dispLimit   = localOpt(S, 'dispLimit',   Inf);
    if isfield(S, 'p0')
        p0 = S.p0;
    else
        p0 = defaultCurveParamInitial(arm);
    end

    ppRoot = fullfile(resultRoot, 'postprocess');
    if ~exist(ppRoot, 'dir'), mkdir(ppRoot); end

    fprintf('Postprocessing curve-parametric result:\n  %s\n', resultRoot);

    [models, analyses] = rebuildCurvePostprocessModels(arm, configs);
    modelRef  = models{1};
    constFull = armConstRingElementIds(modelRef, arm, "full");
    targetVf  = localOpt(curveOpts, 'VolFrac', mean(rhoFullBest));
    rhoMin    = localOpt(curveOpts, 'rhoMin', arm.mma.xminValue);

    rhoFullPipe     = ones(numel(rhoFullBest), 1);
    rhoThreshold05  = double(rhoFullBest > 0.5);
    rhoThreshold05(constFull) = 1.0;
    rhoThreshold05  = max(rhoMin, rhoThreshold05);
    rhoVolumeBinary = makeVolumePreservingBinaryDensity(rhoFullBest, targetVf, constFull);
    rhoVolumeBinary = max(rhoMin, rhoVolumeBinary);

    variants = {
        struct('name', "full_pipe",        'rho', rhoFullPipe);
        struct('name', "threshold_05",     'rho', rhoThreshold05);
        struct('name', "threshold_volume", 'rho', rhoVolumeBinary);
    };

    metricsRows = struct([]);
    for v = 1:numel(variants)
        name = variants{v}.name;
        rho  = variants{v}.rho;
        [metrics, ~] = evaluateLinkedDensityMetrics(analyses, rho, curveOpts);
        printCurveMetrics("Postprocess " + name, metrics, configs);
        metricsRows = [metricsRows; metricRowsForVariant(name, metrics, configs)]; %#ok<AGROW>
    end

    metricsTable = struct2table(metricsRows);
    writetable(metricsTable, fullfile(ppRoot, 'curve_postprocess_metrics.csv'));

    % --- 3-D topology exports (binary variants, no origin colours yet) ------
    fprintf('Exporting topology STL/OBJ (threshold_volume)...\n');
    exportTopology3D(modelRef, rhoVolumeBinary, [], [], ppRoot, 'topology_threshold_volume');
    fprintf('Exporting topology STL/OBJ (threshold_05)...\n');
    exportTopology3D(modelRef, rhoThreshold05,  [], [], ppRoot, 'topology_threshold_05');

    % --- Curve-family origin classification ---------------------------------
    originFull = []; originNames = [];
    if nargoutForBuildCurveLinkedDensity() >= 4
        [~, ~, ~, fields] = buildCurveLinkedDensity(pBest, modelRef, curveOpts);

        % Family fields — single module only.
        exportFamilyFields(modelRef, fields.reference, ppRoot);

        [originRef, originNames, originConfidence] = ...
            classifyCurveFamilyOrigin(fields.reference, 0.85);
        originFull = expandReferenceLabelsToArm(modelRef, originRef);

        % Re-export topology with origin colours.
        fprintf('Exporting coloured topology STL/OBJ...\n');
        exportTopology3D(modelRef, rhoVolumeBinary, originFull, originNames, ...
            ppRoot, 'topology_threshold_volume');
        exportTopology3D(modelRef, rhoThreshold05,  originFull, originNames, ...
            ppRoot, 'topology_threshold_05');

        saveOriginSummary(ppRoot, rhoFullBest, rhoVolumeBinary, ...
            originFull, originNames, originConfidence);
    end

    % --- Initial design single-module export --------------------------------
    fprintf('Exporting initial-design module...\n');
    [~, rhoFullInit] = buildCurveLinkedDensity(p0, modelRef, curveOpts);
    rhoFullInit = rhoFullInit(:);
    rhoVolumeBinaryInit = makeVolumePreservingBinaryDensity(rhoFullInit, targetVf, constFull);
    rhoVolumeBinaryInit = max(rhoMin, rhoVolumeBinaryInit);
    exportTopology3D(modelRef, rhoVolumeBinaryInit, originFull, originNames, ...
        ppRoot, 'initial_topology_module');

    % --- Adversarial config exports -----------------------------------------
    % Re-evaluate with the continuous optimised density (rhoFullBest): that is
    % what the constraint was satisfied with.  rhoVolumeBinary is used only for
    % the 3-D visualisation export.
    if ~isempty(adversarialConfigs)
        adversarialConfigs = reevaluateAdversarialConfigs( ...
            arm, adversarialConfigs, rhoFullBest, curveOpts, stressLimit, dispLimit);
        exportAdversarialConfigTopologies(arm, adversarialConfigs, rhoVolumeBinary, ...
            originFull, originNames, stressLimit, dispLimit, ppRoot);
    end

    writeTextSummary(ppRoot, pBest, metricsTable, targetVf);
    fprintf('Postprocess saved to: %s\n', ppRoot);
end

% =========================================================================
% Local helpers
% =========================================================================

function exportTopology3D(model, rho, originFull, originNames, ppRoot, stem)
% Export full arm: voxel STL + smooth STL + coloured OBJ (if origin provided).
    threshold = 0.5;
    nModules  = [];   % 0 / empty = full arm (all modules)
    hasOrigin = ~isempty(originFull) && ~isempty(originNames);

    exportDensitySTL(model, rho, fullfile(ppRoot, stem + ".stl"), threshold, nModules);

    smoothOpts = struct();
    smoothOpts.nodalSmoothingIters   = 1;
    smoothOpts.surfaceSmoothingIters = 12;
    if hasOrigin
        smoothOpts.labels         = originFull;
        smoothOpts.labelNames     = originNames;
        smoothOpts.coloredObjPath = fullfile(ppRoot, stem + "_colored.obj");
    end
    exportDensitySmoothSTL(model, rho, fullfile(ppRoot, stem + "_smooth.stl"), ...
        threshold, nModules, smoothOpts);

    if hasOrigin
        exportDensityColoredOBJ(model, rho, originFull, originNames, ...
            fullfile(ppRoot, stem + "_voxel_colored.obj"), threshold, nModules);
    end
end

function exportFamilyFields(model, fields, ppRoot)
% Export one-module STL per curve-family field (reference half-segment values only).
    fieldsRoot = fullfile(ppRoot, 'family_fields');
    if ~exist(fieldsRoot, 'dir'), mkdir(fieldsRoot); end
    H          = model.halfSegmentNelems;
    totalElems = model.analysis.getTotalElemsNumber();
    names      = fieldnames(fields);
    for i = 1:numel(names)
        name   = names{i};
        values = fields.(name);
        if numel(values) ~= H, continue; end
        rho = [values(:); zeros(totalElems - H, 1)];
        if max(rho) > 0, rho = rho ./ max(rho); end
        exportDensitySTL(model, rho, fullfile(fieldsRoot, "family_" + name + ".stl"), 0.3, 1);
    end
    fprintf('  Family-field STLs written to: %s\n', fieldsRoot);
end

function adversarialConfigs = reevaluateAdversarialConfigs( ...
        arm, adversarialConfigs, rho, curveOpts, stressLimit, dispLimit)
% Re-evaluate every adversarial config at the final design.
% The violated flag stored at discovery time reflects an intermediate design;
% this corrects it to the actual final-design stress/disp ratios.
    nConfigs = numel(adversarialConfigs);
    fprintf('Re-evaluating %d adversarial configs at final design...\n', nConfigs);

    analyses = cell(nConfigs, 1);
    for ci = 1:nConfigs
        model = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
            arm.res, arm.res_th, arm.alpha, adversarialConfigs(ci).beta, arm.ShapeFn, ...
            false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
        analyses{ci} = model.analysis;
    end

    [metrics, ~] = evaluateLinkedDensityMetrics(analyses, rho, curveOpts);

    for ci = 1:nConfigs
        adversarialConfigs(ci).maxHM       = metrics.maxHM(ci);
        adversarialConfigs(ci).tipUz       = metrics.tipUz(ci);
        adversarialConfigs(ci).stressRatio = metrics.maxHM(ci) / max(stressLimit, eps);
        adversarialConfigs(ci).dispRatio   = abs(metrics.tipUz(ci)) / max(dispLimit, eps);
        adversarialConfigs(ci).violated    = adversarialConfigs(ci).stressRatio > 1.0 || ...
                                             adversarialConfigs(ci).dispRatio   > 1.0;
    end

    nViolated = sum([adversarialConfigs.violated]);
    fprintf('  %d violated, %d satisfied at final design.\n', nViolated, nConfigs - nViolated);

    if nViolated > 0
        fprintf('  Violated configs at final design:\n');
        for ci = 1:nConfigs
            c = adversarialConfigs(ci);
            if ~c.violated, continue; end
            betaStr = sprintf('%6.1f', c.beta);
            fprintf('    [%s]  s_ratio=%.3f', betaStr, c.stressRatio);
            if isfinite(dispLimit)
                fprintf('  d_ratio=%.3f', c.dispRatio);
            end
            fprintf('\n');
        end
    end
end

function exportAdversarialConfigTopologies(arm, adversarialConfigs, rhoVolumeBinary, ...
        originFull, originNames, ~, ~, ppRoot)
    advRoot = fullfile(ppRoot, 'adversarial_configs');
    if ~exist(advRoot, 'dir'), mkdir(advRoot); end
    fprintf('Exporting %d adversarial config modules...\n', numel(adversarialConfigs));
    for ci = 1:numel(adversarialConfigs)
        c = adversarialConfigs(ci);
        model = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
            arm.res, arm.res_th, arm.alpha, c.beta, arm.ShapeFn, ...
            false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
        if c.violated, tag = 'VIOLATED'; else, tag = 'ok'; end
        stem = sprintf('adversarial_%02d_%s', ci, tag);
        exportTopology3D(model, rhoVolumeBinary, originFull, originNames, advRoot, stem);
    end

    rows = struct([]);
    for ci = 1:numel(adversarialConfigs)
        c = adversarialConfigs(ci);
        row.ci          = ci;
        row.stressRatio = c.stressRatio;
        row.dispRatio   = c.dispRatio;
        row.violated    = double(c.violated);
        row.maxHM_MPa   = c.maxHM / 1e6;
        row.tipUz_m     = c.tipUz;
        for ji = 1:numel(c.beta)
            row.(sprintf('beta%d', ji)) = c.beta(ji);
        end
        rows = [rows; row]; %#ok<AGROW>
    end
    if ~isempty(rows)
        writetable(struct2table(rows), fullfile(advRoot, 'adversarial_config_summary.csv'));
    end
    fprintf('Adversarial config exports saved to: %s\n', advRoot);
end

function [models, analyses] = rebuildCurvePostprocessModels(arm, configs)
    nConfigs = numel(configs);
    models   = cell(nConfigs, 1);
    analyses = cell(nConfigs, 1);
    for k = 1:nConfigs
        cfg = configs{k};
        model = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
            arm.res, arm.res_th, arm.alpha, cfg.betas, arm.ShapeFn, false, ...
            arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
        if k > 1
            assertLinkedArmLayoutCompatible(model, models{1}, cfg.name);
        end
        models{k}   = model;
        analyses{k} = model.analysis;
    end
end

function rows = metricRowsForVariant(variantName, metrics, configs)
    rows = repmat(struct(), numel(configs), 1);
    for k = 1:numel(configs)
        rows(k).variant         = string(variantName);
        rows(k).configName      = string(configs{k}.name);
        rows(k).configLabel     = string(configs{k}.label);
        rows(k).volumeFraction  = metrics.volumeFraction;
        rows(k).maxHM           = metrics.maxHM(k);
        rows(k).stressAggregate = metrics.stressAggregateByConfig(k);
        rows(k).tipUz           = metrics.tipUz(k);
        rows(k).maxDisp         = metrics.maxDisp(k);
    end
end

function saveOriginSummary(resultRoot, rhoContinuous, rhoBinary, originFull, originNames, originConfidence)
    names = ["mixed"; originNames(:)];
    rows  = repmat(struct(), numel(names), 1);
    selectedContinuous = rhoContinuous(:) > 0.5;
    selectedBinary     = rhoBinary(:) > 0.5;
    originRefRepeatedConfidence = expandConfidence(originConfidence, numel(originFull));
    for i = 1:numel(names)
        label = i - 1;
        rows(i).origin             = names(i);
        rows(i).continuousCount    = nnz(selectedContinuous & originFull(:) == label);
        rows(i).binaryCount        = nnz(selectedBinary     & originFull(:) == label);
        rows(i).continuousFraction = rows(i).continuousCount / max(1, nnz(selectedContinuous));
        rows(i).binaryFraction     = rows(i).binaryCount     / max(1, nnz(selectedBinary));
        mask = selectedBinary & originFull(:) == label;
        rows(i).meanConfidence     = mean(originRefRepeatedConfidence(mask), 'omitnan');
    end
    writetable(struct2table(rows), fullfile(resultRoot, 'curve_origin_summary.csv'));
end

function confidenceFull = expandConfidence(confidenceRef, nFull)
    H = numel(confidenceRef);
    nCopies = nFull / H;
    confidenceFull = confidenceRef(:);
    for k = 1:(nCopies/2 - 1)
        confidenceFull = [confidenceFull; flip(confidenceRef(:)); confidenceRef(:)]; %#ok<AGROW>
    end
    confidenceFull = [confidenceFull; flip(confidenceRef(:))];
end

function writeTextSummary(resultRoot, params, metricsTable, targetVf)
    fid     = fopen(fullfile(resultRoot, 'curve_postprocess_summary.txt'), 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'Curve-parametric postprocess summary\n');
    fprintf(fid, 'Target volume fraction: %.6f\n\n', targetVf);
    fprintf(fid, 'Best parameters:\n');
    names = curveParamNames();
    for i = 1:numel(names)
        fprintf(fid, '  %-20s %.8g\n', names{i}, params.(names{i}));
    end
    fprintf(fid, '\nMetrics by variant/config:\n');
    for i = 1:height(metricsTable)
        fprintf(fid, '  %-18s %-12s maxHM=% .6e stressAgg=% .6e tipUz=% .6e maxDisp=% .6e vf=%.5f\n', ...
            char(metricsTable.variant(i)), char(metricsTable.configName(i)), ...
            metricsTable.maxHM(i), metricsTable.stressAggregate(i), ...
            metricsTable.tipUz(i), metricsTable.maxDisp(i), metricsTable.volumeFraction(i));
    end
end

function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name)
        value = s.(name);
    else
        value = defaultValue;
    end
end

function n = nargoutForBuildCurveLinkedDensity()
    n = nargout('buildCurveLinkedDensity');
end

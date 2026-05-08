function postprocessCurveParamResult(resultRoot)
% postprocessCurveParamResult  Reanalyse and plot curve-parametric designs.
%
% Creates thresholded topology plots, an unwrapped reference half-segment
% topology map, binary reanalysis metrics (full-pipe reference, optimized
% threshold_05, threshold_volume), and comparison charts.

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

    pBest = S.pBest;
    rhoFullBest = S.rhoFullBest(:);
    curveOpts = S.curveOpts;
    arm = S.arm;
    if isRobustResult
        % Use a single nominal (straight) config as reference for plots.
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
    modelRef = models{1};
    H = modelRef.halfSegmentNelems;
    constFull = armConstRingElementIds(modelRef, arm, "full");
    constRef = armConstRingElementIds(modelRef, arm, "linked");
    targetVf = localOpt(curveOpts, 'VolFrac', mean(rhoFullBest));
    rhoMin = localOpt(curveOpts, 'rhoMin', arm.mma.xminValue);

    rhoFullPipe = ones(numel(rhoFullBest), 1);
    rhoThreshold05 = double(rhoFullBest > 0.5);
    rhoThreshold05(constFull) = 1.0;
    rhoThreshold05 = max(rhoMin, rhoThreshold05);
    rhoVolumeBinary = makeVolumePreservingBinaryDensity(rhoFullBest, targetVf, constFull);
    rhoVolumeBinary = max(rhoMin, rhoVolumeBinary);

    variants = {
        struct('name', "full_pipe",       'rho', rhoFullPipe);
        struct('name', "threshold_05",    'rho', rhoThreshold05);
        struct('name', "threshold_volume",'rho', rhoVolumeBinary);
    };

    metricsRows = struct([]);
    metricsVol = [];
    for v = 1:numel(variants)
        name = variants{v}.name;
        rho = variants{v}.rho;
        [metrics, ~] = evaluateLinkedDensityMetrics(analyses, rho, curveOpts);
        printCurveMetrics("Postprocess " + name, metrics, configs);
        metricsRows = [metricsRows; metricRowsForVariant(name, metrics, configs)]; %#ok<AGROW>
        if name == "threshold_volume"
            metricsVol = metrics;
        end
    end

    metricsTable = struct2table(metricsRows);
    writetable(metricsTable, fullfile(ppRoot, 'curve_postprocess_metrics.csv'));

    rhoRefContinuous = rhoFullBest(1:H);
    rhoRefVolumeBinary = rhoVolumeBinary(1:H);

    saveTopologyFigure(plotCurveTopology(modelRef, rhoVolumeBinary, ...
        sprintf('Volume-preserving threshold, vf=%.3f', mean(rhoVolumeBinary > rhoMin)), "threshold"), ...
        ppRoot, 'rho_full_threshold_volume');
    saveTopologyFigure(plotCurveTopology(modelRef, rhoVolumeBinary, ...
        sprintf('Real topology, volume threshold, vf=%.3f', mean(rhoVolumeBinary > rhoMin)), "real"), ...
        ppRoot, 'topology_full_real_threshold_volume');
    exportCurveTopologyPatch(modelRef, rhoVolumeBinary, ...
        fullfile(ppRoot, 'topology_full_real_threshold_volume_patch.mat'), 0.5);

    saveTopologyFigure(plotCurveTopology(modelRef, rhoThreshold05, ...
        sprintf('Real topology, rho > 0.5, vf=%.3f', mean(rhoThreshold05 > rhoMin)), "real"), ...
        ppRoot, 'topology_full_real_threshold_05');
    exportCurveTopologyPatch(modelRef, rhoThreshold05, ...
        fullfile(ppRoot, 'topology_full_real_threshold_05_patch.mat'), 0.5);

    saveTopologyFigure(plotCurveTopology(modelRef, [rhoRefVolumeBinary; zeros(numel(rhoFullBest)-H, 1)], ...
        'Reference half-segment thresholded topology', "threshold"), ...
        ppRoot, 'rho_reference_threshold_volume');
    saveTopologyFigure(plotCurveTopology(modelRef, [rhoRefVolumeBinary; zeros(numel(rhoFullBest)-H, 1)], ...
        'Reference half-segment real topology', "real"), ...
        ppRoot, 'topology_reference_real_threshold_volume');

    saveTopologyFigure(plotCurveUnwrappedTopology(modelRef, rhoRefContinuous, pBest, ...
        'Unwrapped reference density with curve-family overlays'), ...
        ppRoot, 'topology_unwrapped_reference');
    saveTopologyFigure(plotCurveUnwrappedTopology(modelRef, rhoRefVolumeBinary, pBest, ...
        'Unwrapped reference thresholded topology'), ...
        ppRoot, 'topology_unwrapped_reference_threshold');

    fprintf('Generating initial-parameter comparison figures...\n');
    [rhoRefInit, rhoFullInit] = buildCurveLinkedDensity(p0, modelRef, curveOpts);
    rhoFullInit = rhoFullInit(:);
    rhoVolumeBinaryInit = makeVolumePreservingBinaryDensity(rhoFullInit, targetVf, constFull);
    rhoVolumeBinaryInit = max(rhoMin, rhoVolumeBinaryInit);
    rhoRefVolumeBinaryInit = rhoVolumeBinaryInit(1:H);

    saveTopologyFigure(plotCurveTopology(modelRef, rhoVolumeBinaryInit, ...
        sprintf('Initial: real topology, volume threshold, vf=%.3f', mean(rhoVolumeBinaryInit > rhoMin)), "real"), ...
        ppRoot, 'initial_topology_full_real_threshold_volume');
    saveTopologyFigure(plotCurveTopology(modelRef, [rhoRefVolumeBinaryInit; zeros(numel(rhoFullInit)-H, 1)], ...
        'Initial: reference half-segment real topology', "real"), ...
        ppRoot, 'initial_topology_reference_real_threshold_volume');
    saveTopologyFigure(plotCurveUnwrappedTopology(modelRef, rhoRefInit(:), p0, ...
        'Initial: unwrapped reference density with curve-family overlays'), ...
        ppRoot, 'initial_topology_unwrapped_reference');
    saveTopologyFigure(plotCurveUnwrappedTopology(modelRef, rhoRefVolumeBinaryInit, p0, ...
        'Initial: unwrapped reference thresholded topology'), ...
        ppRoot, 'initial_topology_unwrapped_reference_threshold');

    plotMetricSet(metricsTable, ppRoot);

    if nargoutForBuildCurveLinkedDensity() >= 4
        [~, ~, ~, fields] = buildCurveLinkedDensity(pBest, modelRef, curveOpts);
        plotFamilyFields(modelRef, fields.reference, ppRoot);
        [originRef, originNames, originConfidence] = classifyCurveFamilyOrigin(fields.reference, 0.85);
        originFull = expandReferenceLabelsToArm(modelRef, originRef);

        saveTopologyFigure(plotCurveTopologyByOrigin(modelRef, rhoVolumeBinary, originFull, originNames, ...
            'Real topology colored by dominant curve family: volume threshold', 0.5), ...
            ppRoot, 'topology_full_real_origin_threshold_volume');
        saveTopologyFigure(plotCurveTopologyByOrigin(modelRef, rhoThreshold05, originFull, originNames, ...
            'Real topology colored by dominant curve family: rho > 0.5', 0.5), ...
            ppRoot, 'topology_full_real_origin_threshold_05');

        % Per-configuration topology plots (all arm poses).
        for k = 1:numel(models)
            cfgLabel = string(configs{k}.label);
            cfgName  = string(configs{k}.name);
            saveTopologyFigure(plotCurveTopologyByOrigin(models{k}, rhoVolumeBinary, originFull, originNames, ...
                'Real topology: ' + cfgLabel, 0.5), ...
                ppRoot, 'topology_origin_' + cfgName);
        end

        % Per-configuration origin topology — one image per pose.
        for k = 1:numel(models)
            cfgLabel = string(configs{k}.label);
            cfgName  = string(configs{k}.name);
            titleLines = {char(cfgLabel)};
            if ~isempty(metricsVol)
                constraintStr = '';
                if isfinite(stressLimit)
                    sVal = metricsVol.maxHM(k);
                    if sVal <= stressLimit, sTag = 'OK'; else, sTag = 'FAIL'; end
                    constraintStr = sprintf('sigma: %.2f / %.2f MPa [%s]', sVal/1e6, stressLimit/1e6, sTag);
                end
                if isfinite(dispLimit)
                    uVal = abs(metricsVol.tipUz(k));
                    if uVal <= dispLimit, uTag = 'OK'; else, uTag = 'FAIL'; end
                    uStr = sprintf('u: %.1f / %.1f mm [%s]', uVal*1e3, dispLimit*1e3, uTag);
                    if isempty(constraintStr)
                        constraintStr = uStr;
                    else
                        constraintStr = [constraintStr, '   ', uStr];
                    end
                end
                if ~isempty(constraintStr)
                    titleLines{2} = constraintStr;
                end
            end
            saveTopologyFigure(plotCurveTopologyByOrigin(models{k}, rhoVolumeBinary, originFull, originNames, ...
                titleLines, 0.5), ...
                ppRoot, 'topology_all_configs_origin_' + cfgName);
        end
        saveTopologyFigure(plotCurveUnwrappedBarsByOrigin(modelRef, rhoVolumeBinary(1:H), ...
            originRef, originNames, pBest, ...
            'Unwrapped bar topology colored by dominant curve family', 0.5), ...
            ppRoot, 'topology_unwrapped_bars_origin_threshold_volume');
        saveOriginSummary(ppRoot, rhoFullBest, rhoVolumeBinary, originFull, originNames, originConfidence);

        % ---- Adversarial config topology plots (robust results only). ----
        if ~isempty(adversarialConfigs)
            plotAdversarialConfigTopologies(arm, adversarialConfigs, rhoVolumeBinary, ...
                originFull, originNames, stressLimit, dispLimit, ppRoot);
        end
    end

    writeTextSummary(ppRoot, pBest, metricsTable, targetVf);
    fprintf('Postprocess saved to: %s\n', ppRoot);
end

function [models, analyses] = rebuildCurvePostprocessModels(arm, configs)
    nConfigs = numel(configs);
    models = cell(nConfigs, 1);
    analyses = cell(nConfigs, 1);
    for k = 1:nConfigs
        cfg = configs{k};
        model = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
            arm.res, arm.res_th, arm.alpha, cfg.betas, arm.ShapeFn, false, ...
            arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
        if k > 1
            assertLinkedArmLayoutCompatible(model, models{1}, cfg.name);
        end
        models{k} = model;
        analyses{k} = model.analysis;
    end
end

function rows = metricRowsForVariant(variantName, metrics, configs)
    rows = repmat(struct(), numel(configs), 1);
    for k = 1:numel(configs)
        rows(k).variant = string(variantName);
        rows(k).configName = string(configs{k}.name);
        rows(k).configLabel = string(configs{k}.label);
        rows(k).volumeFraction = metrics.volumeFraction;
        rows(k).maxHM = metrics.maxHM(k);
        rows(k).stressAggregate = metrics.stressAggregateByConfig(k);
        rows(k).tipUz = metrics.tipUz(k);
        rows(k).maxDisp = metrics.maxDisp(k);
    end
end

function saveTopologyFigure(fig, resultRoot, stem)
    set(fig, 'Visible', 'on');
    exportgraphics(fig, fullfile(resultRoot, stem + ".png"), 'Resolution', 220);
    savefig(fig, fullfile(resultRoot, stem + ".fig"));
end

function plotMetricSet(metricsTable, resultRoot)
    metrics = {'maxHM', 'stressAggregate', 'tipUz', 'maxDisp'};
    for i = 1:numel(metrics)
        metric = metrics{i};
        fig = plotCurveMetricBars(metricsTable, metric, "Curve postprocess: " + metric);
        saveTopologyFigure(fig, resultRoot, "bar_" + metric);
    end
end

function plotFamilyFields(model, fields, resultRoot)
    names = fieldnames(fields);
    H = model.halfSegmentNelems;
    for i = 1:numel(names)
        name = names{i};
        values = fields.(name);
        if numel(values) ~= H
            continue;
        end
        rho = [values(:); zeros(model.analysis.getTotalElemsNumber() - H, 1)];
        if max(rho) > 0
            rho = rho ./ max(rho);
        end
        fig = plotCurveTopology(model, rho, "Curve family field: " + string(name), "continuous");
        saveTopologyFigure(fig, resultRoot, "family_" + string(name));
    end
end

function saveOriginSummary(resultRoot, rhoContinuous, rhoBinary, originFull, originNames, originConfidence)
    names = ["mixed"; originNames(:)];
    rows = repmat(struct(), numel(names), 1);
    selectedContinuous = rhoContinuous(:) > 0.5;
    selectedBinary = rhoBinary(:) > 0.5;
    originRefRepeatedConfidence = expandConfidence(originConfidence, numel(originFull));

    for i = 1:numel(names)
        label = i - 1;
        rows(i).origin = names(i);
        rows(i).continuousCount = nnz(selectedContinuous & originFull(:) == label);
        rows(i).binaryCount = nnz(selectedBinary & originFull(:) == label);
        rows(i).continuousFraction = rows(i).continuousCount / max(1, nnz(selectedContinuous));
        rows(i).binaryFraction = rows(i).binaryCount / max(1, nnz(selectedBinary));
        mask = selectedBinary & originFull(:) == label;
        rows(i).meanConfidence = mean(originRefRepeatedConfidence(mask), 'omitnan');
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
    fid = fopen(fullfile(resultRoot, 'curve_postprocess_summary.txt'), 'w');
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
            char(metricsTable.variant(i)), char(metricsTable.configName(i)), metricsTable.maxHM(i), ...
            metricsTable.stressAggregate(i), metricsTable.tipUz(i), metricsTable.maxDisp(i), ...
            metricsTable.volumeFraction(i));
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

function plotAdversarialConfigTopologies(arm, adversarialConfigs, rhoVolumeBinary, ...
    originFull, originNames, stressLimit, dispLimit, ppRoot)
% Plot one topology image per adversarial candidate (visually distinct titles).
% Violated candidates are labelled [VIOLATED]; non-violated are [OK].

    advRoot = fullfile(ppRoot, 'adversarial_configs');
    if ~exist(advRoot, 'dir'), mkdir(advRoot); end

    fprintf('Plotting %d adversarial config topologies...\n', numel(adversarialConfigs));
    for ci = 1:numel(adversarialConfigs)
        c = adversarialConfigs(ci);
        betaVec = c.beta;

        % Build model for this pose.
        model = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
            arm.res, arm.res_th, arm.alpha, betaVec, arm.ShapeFn, ...
            false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);

        % Build two-line title.
        if c.violated, statusTag = 'VIOLATED'; else, statusTag = 'ok'; end
        line1 = sprintf('Adversarial config %d [%s]', ci, statusTag);

        constraintStr = '';
        if isfinite(stressLimit)
            if c.stressRatio > 1, sTag = 'FAIL'; else, sTag = 'OK'; end
            constraintStr = sprintf('sigma: %.2f MPa (ratio %.2f) [%s]', ...
                c.maxHM/1e6, c.stressRatio, sTag);
        end
        if isfinite(dispLimit)
            if c.dispRatio > 1, uTag = 'FAIL'; else, uTag = 'OK'; end
            uStr = sprintf('u: %.2f mm (ratio %.2f) [%s]', ...
                abs(c.tipUz)*1e3, c.dispRatio, uTag);
            if isempty(constraintStr)
                constraintStr = uStr;
            else
                constraintStr = [constraintStr, '   ', uStr];
            end
        end
        betaStr = sprintf('beta=[%s]', num2str(betaVec, '%g '));

        titleLines = {line1, constraintStr, betaStr};
        titleLines = titleLines(~cellfun('isempty', titleLines));

        stem = sprintf('adversarial_%02d_%s', ci, statusTag);
        saveTopologyFigure(plotCurveTopologyByOrigin(model, rhoVolumeBinary, ...
            originFull, originNames, titleLines, 0.5), advRoot, stem);
    end

    % Write summary CSV for adversarial configs.
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
        writetable(struct2table(rows), fullfile(advRoot, 'adversarial_config_summary.csv'));
    end
    fprintf('Adversarial config plots saved to: %s\n', advRoot);
end

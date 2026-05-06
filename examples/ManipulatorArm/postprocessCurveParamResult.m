function postprocessCurveParamResult(resultRoot)
% postprocessCurveParamResult  Reanalyse and plot curve-parametric designs.
%
% Creates continuous and thresholded topology plots, an unwrapped reference
% half-segment topology map, binary reanalysis metrics, and comparison charts.

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
    required = {'pBest', 'rhoFullBest', 'configs', 'curveOpts', 'arm'};
    for i = 1:numel(required)
        assert(isfield(S, required{i}), 'result.mat lacks "%s".', required{i});
    end

    pBest = S.pBest;
    rhoFullBest = S.rhoFullBest(:);
    configs = S.configs;
    curveOpts = S.curveOpts;
    arm = S.arm;

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

    rhoContinuous = rhoFullBest;
    rhoThreshold05 = double(rhoFullBest > 0.5);
    rhoThreshold05(constFull) = 1.0;
    rhoThreshold05 = max(rhoMin, rhoThreshold05);
    rhoVolumeBinary = makeVolumePreservingBinaryDensity(rhoFullBest, targetVf, constFull);
    rhoVolumeBinary = max(rhoMin, rhoVolumeBinary);

    variants = {
        struct('name', "continuous", 'rho', rhoContinuous);
        struct('name', "threshold_05", 'rho', rhoThreshold05);
        struct('name', "threshold_volume", 'rho', rhoVolumeBinary);
    };

    metricsRows = struct([]);
    for v = 1:numel(variants)
        name = variants{v}.name;
        rho = variants{v}.rho;
        [metrics, ~] = evaluateLinkedDensityMetrics(analyses, rho, curveOpts);
        printCurveMetrics("Postprocess " + name, metrics, configs);
        metricsRows = [metricsRows; metricRowsForVariant(name, metrics, configs)]; %#ok<AGROW>
    end

    metricsTable = struct2table(metricsRows);
    writetable(metricsTable, fullfile(ppRoot, 'curve_postprocess_metrics.csv'));

    rhoRefContinuous = rhoContinuous(1:H);
    rhoRefVolumeBinary = rhoVolumeBinary(1:H);

    saveTopologyFigure(plotCurveTopology(modelRef, rhoContinuous, ...
        sprintf('Continuous curve density, vf=%.3f', mean(rhoContinuous)), "continuous"), ...
        ppRoot, 'rho_full_continuous');
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

    saveTopologyFigure(plotCurveTopology(modelRef, [rhoRefContinuous; zeros(numel(rhoContinuous)-H, 1)], ...
        'Reference half-segment continuous density', "continuous"), ...
        ppRoot, 'rho_reference_continuous');
    saveTopologyFigure(plotCurveTopology(modelRef, [rhoRefVolumeBinary; zeros(numel(rhoVolumeBinary)-H, 1)], ...
        'Reference half-segment thresholded topology', "threshold"), ...
        ppRoot, 'rho_reference_threshold_volume');
    saveTopologyFigure(plotCurveTopology(modelRef, [rhoRefVolumeBinary; zeros(numel(rhoVolumeBinary)-H, 1)], ...
        'Reference half-segment real topology', "real"), ...
        ppRoot, 'topology_reference_real_threshold_volume');

    saveTopologyFigure(plotCurveUnwrappedTopology(modelRef, rhoRefContinuous, pBest, ...
        'Unwrapped reference density with curve-family overlays'), ...
        ppRoot, 'topology_unwrapped_reference');
    saveTopologyFigure(plotCurveUnwrappedTopology(modelRef, rhoRefVolumeBinary, pBest, ...
        'Unwrapped reference thresholded topology'), ...
        ppRoot, 'topology_unwrapped_reference_threshold');

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
        saveTopologyFigure(plotCurveUnwrappedBarsByOrigin(modelRef, rhoVolumeBinary(1:H), ...
            originRef, originNames, pBest, ...
            'Unwrapped bar topology colored by dominant curve family', 0.5), ...
            ppRoot, 'topology_unwrapped_bars_origin_threshold_volume');
        saveOriginSummary(ppRoot, rhoFullBest, rhoVolumeBinary, originFull, originNames, originConfidence);
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

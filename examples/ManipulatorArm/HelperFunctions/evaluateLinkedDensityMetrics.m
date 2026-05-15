function [metrics, raw] = evaluateLinkedDensityMetrics(analyses, rhoFull, opts)
% evaluateLinkedDensityMetrics  Smooth stress and vertical displacement metrics.

    nConfigs = numel(analyses);
    rhoFull = rhoFull(:);
    penal = localOpt(opts, 'penal', 3.0);
    pElem = localOpt(opts, 'pElem', 12.0);
    pConfigStress = localOpt(opts, 'pConfigStress', 4.0);
    pConfigDisp = localOpt(opts, 'pConfigDisp', 4.0);
    objectiveMode = char(string(localOpt(opts, 'objectiveMode', "stress")));

    maxHM = zeros(nConfigs, 1);
    stressAgg = zeros(nConfigs, 1);
    tipUz = zeros(nConfigs, 1);
    maxDisp = zeros(nConfigs, 1);
    stressValues = cell(nConfigs, 1);

    xPenal = rhoFull .^ penal;
    parfor k = 1:nConfigs
        analysis = analyses{k};
        analysis.solveWeighted(xPenal);
        analysis.computeElementResults(xPenal);

        fe = analysis.felems{1};
        hmIdx = find(fe.results.names == "sHM", 1);
        hm = squeeze(fe.results.gp.all(hmIdx, :, :));
        hm = hm(:);
        hm = hm(isfinite(hm));
        if isempty(hm)
            hm = 0;
        end

        maxHM(k) = max(hm);
        stressAgg(k) = smoothPnorm(hm, pElem);
        stressValues{k} = hm;

        uzIdx = analysis.findDOFsIndices("uz");
        loadedNodes = find(any(abs(analysis.Pnodal) > 0, 2));
        if isempty(loadedNodes)
            loadedNodes = size(analysis.qnodal, 1);
        end
        tipUz(k) = mean(analysis.qnodal(loadedNodes, uzIdx));

        dispIdx = analysis.findDOFsIndices(["ux" "uy" "uz"]);
        maxDisp(k) = max(vecnorm(analysis.qnodal(:, dispIdx), 2, 2));
    end

    stressRef = localOpt(opts, 'stressRef', ones(nConfigs, 1));
    dispRef = localOpt(opts, 'dispRef', ones(nConfigs, 1));
    stressRef = max(abs(stressRef(:)), eps);
    dispRef = max(abs(dispRef(:)), eps);
    configWeights = localOpt(opts, 'configWeights', ones(nConfigs, 1) / nConfigs);
    configWeights = configWeights(:);
    assert(numel(configWeights) == nConfigs, ...
        'configWeights must have one entry per configuration.');
    configWeights = max(0, configWeights);
    if sum(configWeights) <= 0
        configWeights = ones(nConfigs, 1) / nConfigs;
    else
        configWeights = configWeights / sum(configWeights);
    end

    stressObjective = weightedSmoothPnorm(stressAgg ./ stressRef, configWeights, pConfigStress);
    dispObjective = weightedSmoothPnorm(abs(tipUz) ./ dispRef, configWeights, pConfigDisp);

    switch objectiveMode
        case 'stress'
            objective = stressObjective;
        case {'uz', 'disp', 'displacement'}
            objective = dispObjective;
        case 'combined'
            stressWeight = localOpt(opts, 'stressWeight', 0.5);
            dispWeight = localOpt(opts, 'dispWeight', 0.5);
            objective = stressWeight * stressObjective + dispWeight * dispObjective;
        case 'volume'
            objective = mean(rhoFull);
        otherwise
            error('Unknown objectiveMode "%s".', objectiveMode);
    end

    metrics = struct();
    metrics.objective = objective;
    metrics.stressObjective = stressObjective;
    metrics.dispObjective = dispObjective;
    metrics.maxHM = maxHM;
    metrics.stressAggregateByConfig = stressAgg;
    metrics.tipUz = tipUz;
    metrics.maxDisp = maxDisp;
    metrics.volumeFraction = mean(rhoFull);

    raw = struct();
    raw.stressValues = stressValues;
end

function y = smoothPnorm(x, p)
    x = abs(x(:));
    if isempty(x)
        y = 0;
        return;
    end
    m = max(x);
    if m <= 0
        y = 0;
    else
        y = m * mean((x / m) .^ p) .^ (1 / p);
    end
end

function y = weightedSmoothPnorm(x, w, p)
    x = abs(x(:));
    w = w(:);
    if isempty(x)
        y = 0;
        return;
    end
    m = max(x);
    if m <= 0
        y = 0;
    else
        y = m * sum(w .* (x / m) .^ p) .^ (1 / p);
    end
end

function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name)
        value = s.(name);
    else
        value = defaultValue;
    end
end

function [S, maxStress] = evaluateStressSetOnly( ...
        analyses, x, penal, stressPNorm, q, Starget, useParallel, nStressClusters)
% EVALUATESTRESSSETONLY  Forward-only stress aggregate evaluation (no gradient).
%
%   [S, maxStress] = evaluateStressSetOnly(analyses, x, penal, stressPNorm, q,
%                                          Starget, useParallel, nStressClusters)
%
%   Replicates the forward pass of evaluateStressSetLocal /
%   computeStressAggregateAndGradient from solveSIMPVolumeStressMMA without
%   the adjoint gradient.  Used for post-processing topology candidates.
%
%   Inputs:
%     analyses         {nConfigs x 1}  FEAnalysis objects
%     x                [nElems x 1]  physical density (full-arm)
%     penal            SIMP penalisation exponent
%     stressPNorm      p-norm aggregation exponent
%     q                qp-relaxation exponent
%     Starget          [nConfigs*nStressClusters x 1] stress targets; [] for unit
%     useParallel      logical
%     nStressClusters  number of stress-level clusters per config
%
%   Outputs:
%     S          [nConfigs*nStressClusters x 1]  p-norm stress aggregates
%     maxStress  [nConfigs*nStressClusters x 1]  max Huber-Mises per cluster

    nConfigs     = numel(analyses);
    nConstraints = nConfigs * nStressClusters;
    S            = zeros(nConstraints, 1);
    maxStress    = zeros(nConstraints, 1);

    useParallel = useParallel && license('test', 'Distrib_Computing_Toolbox');
    if useParallel
        SCell         = cell(nConfigs, 1);
        maxStressCell = cell(nConfigs, 1);
        parfor k = 1:nConfigs
            idx     = (k - 1) * nStressClusters + (1:nStressClusters);
            targetK = [];
            if ~isempty(Starget), targetK = Starget(idx); end
            [SCell{k}, maxStressCell{k}] = computeStressForward( ...
                analyses{k}, x, penal, stressPNorm, q, targetK, nStressClusters);
        end
        for k = 1:nConfigs
            idx           = (k - 1) * nStressClusters + (1:nStressClusters);
            S(idx)         = SCell{k};
            maxStress(idx) = maxStressCell{k};
        end
    else
        for k = 1:nConfigs
            idx     = (k - 1) * nStressClusters + (1:nStressClusters);
            targetK = [];
            if ~isempty(Starget), targetK = Starget(idx); end
            [S(idx), maxStress(idx)] = computeStressForward( ...
                analyses{k}, x, penal, stressPNorm, q, targetK, nStressClusters);
        end
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function [S, maxStress] = computeStressForward( ...
        analysis, x, penal, stressPNorm, q, Starget, nStressClusters)
    nElems = analysis.getTotalElemsNumber();
    x      = x(:);
    xPenal = x .^ penal;
    analysis.solveWeighted(xPenal, false);
    analysis.computeElementResults(xPenal);

    sigma   = elementGPHuberMisesStress(analysis, nElems);
    relaxed = (max(x, eps) .^ q) .* sigma;

    if isempty(Starget)
        target = ones(nStressClusters, 1);
    else
        target = max(Starget(:), eps);
    end

    clusters  = stressLevelClusters(relaxed, nStressClusters);
    S         = zeros(nStressClusters, 1);
    maxStress = zeros(nStressClusters, 1);

    for c = 1:nStressClusters
        elemIds      = clusters{c};
        maxStress(c) = max(sigma(elemIds));
        ratio        = relaxed(elemIds) / target(c);
        meanPower    = mean(ratio .^ stressPNorm);
        S(c)         = target(c) * meanPower ^ (1.0 / stressPNorm);
    end
end

% -------------------------------------------------------------------------
function sigma = elementGPHuberMisesStress(analysis, nElems)
    sigma       = zeros(nElems, 1);
    elemIndices = analysis.getElemIndices();
    for i = 1:numel(analysis.felems)
        fe = analysis.felems{i};
        if ~isfield(fe.results, 'gp') || ~isfield(fe.results.gp, 'stress')
            continue;
        end
        s  = fe.results.gp.stress;   % (nElemsI, nip, 6)
        s1 = s(:,:,1); s2 = s(:,:,2); s3 = s(:,:,3);
        s4 = s(:,:,4); s5 = s(:,:,5); s6 = s(:,:,6);
        hmGP    = sqrt(0.5*((s1-s2).^2 + (s2-s3).^2 + (s3-s1).^2) + ...
                       3*(s4.^2+s5.^2+s6.^2));
        elemIds = elemIndices{i};
        sigma(elemIds) = mean(hmGP, 2);
    end
end

% -------------------------------------------------------------------------
function clusters = stressLevelClusters(stressMeasure, nStressClusters)
    [~, order]      = sort(stressMeasure(:), 'descend');
    nElems          = numel(order);
    nStressClusters = min(max(1, nStressClusters), nElems);
    clusters        = cell(nStressClusters, 1);
    edges           = round(linspace(0, nElems, nStressClusters + 1));
    for c = 1:nStressClusters
        ids = order((edges(c) + 1):edges(c + 1));
        if isempty(ids), ids = order(1); end
        clusters{c} = ids;
    end
end

function [J, gradJ, C0, C] = evaluateObjectiveAndGradient(analyses, x, penal, pAgg, weights, C0, useParallel)
    if nargin < 7
        useParallel = false;
    end
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

    if isempty(C0)
        C0 = max(C, eps);
    end

    Cagg = C ./ C0;
    weightedSum = sum(weights(:) .* (Cagg .^ pAgg));
    J = weightedSum ^ (1.0 / pAgg);

    gradJ = zeros(numel(x), 1);
    if J > eps
        for k = 1:nConfigs
            gradJ = gradJ + weights(k) * Cagg(k)^(pAgg - 1) * dC(:, k) / C0(k);
        end
        gradJ = J^(1 - pAgg) * gradJ;
    end
end

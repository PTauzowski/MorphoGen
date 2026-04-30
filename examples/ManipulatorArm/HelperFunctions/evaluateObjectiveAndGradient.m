function [J, gradJ, C0, C] = evaluateObjectiveAndGradient(analyses, x, penal, pAgg, weights, C0)
    nConfigs = numel(analyses);
    C = zeros(nConfigs, 1);
    dC = zeros(numel(x), nConfigs);

    for k = 1:nConfigs
        [C(k), dC(:, k)] = computeComplianceAndGradient(analyses{k}, x, penal);
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

function [J, C] = evaluateObjectiveOnly(analyses, x, penal, pAgg, weights, C0)
    nConfigs = numel(analyses);
    C = zeros(nConfigs, 1);

    for k = 1:nConfigs
        C(k) = computeComplianceOnly(analyses{k}, x, penal);
    end

    assert(~isempty(C0), 'Compliance normalization C0 must be provided.');
    Cagg = C ./ C0;
    weightedSum = sum(weights(:) .* (Cagg .^ pAgg));
    J = weightedSum ^ (1.0 / pAgg);
end

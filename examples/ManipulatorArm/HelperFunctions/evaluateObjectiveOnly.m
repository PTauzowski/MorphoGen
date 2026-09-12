function [J, C] = evaluateObjectiveOnly(analyses, x, penal, pAgg, weights, C0, useParallel)
    if nargin < 7
        useParallel = false;
    end
    nConfigs = numel(analyses);
    C = zeros(nConfigs, 1);

    useParallel = useParallel && license('test', 'Distrib_Computing_Toolbox');
    if useParallel
        parfor k = 1:nConfigs
            C(k) = computeComplianceOnly(analyses{k}, x, penal);
        end
    else
        for k = 1:nConfigs
            C(k) = computeComplianceOnly(analyses{k}, x, penal);
        end
    end

    assert(~isempty(C0), 'Compliance normalization C0 must be provided.');
    Cagg = C ./ C0;
    weightedSum = sum(weights(:) .* (Cagg .^ pAgg));
    J = weightedSum ^ (1.0 / pAgg);
end

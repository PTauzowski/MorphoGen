function u = curveParamsToUnconstrained(params, bounds)
% curveParamsToUnconstrained  Logistic inverse transform for bounded params.

    names = curveParamNames();
    u = zeros(numel(names), 1);
    epsClamp = 1.0e-8;
    for i = 1:numel(names)
        name = names{i};
        b = bounds.(name);
        lo = b(1);
        hi = b(2);
        t = (params.(name) - lo) / (hi - lo);
        t = min(1 - epsClamp, max(epsClamp, t));
        u(i) = log(t / (1 - t));
    end
end

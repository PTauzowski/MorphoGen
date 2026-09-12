function params = unconstrainedToCurveParams(u, bounds)
% unconstrainedToCurveParams  Logistic transform from R^n to bounded params.

    names = curveParamNames();
    assert(numel(u) == numel(names), 'Expected %d curve parameters.', numel(names));
    params = struct();
    for i = 1:numel(names)
        name = names{i};
        b = bounds.(name);
        lo = b(1);
        hi = b(2);
        t = 1 / (1 + exp(-u(i)));
        params.(name) = lo + (hi - lo) * t;
    end
end

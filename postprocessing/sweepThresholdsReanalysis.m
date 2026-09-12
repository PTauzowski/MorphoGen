function result = sweepThresholdsReanalysis(x, analyses, thresholds, opts)
% SWEEPTHRESHOLDSREANALYSIS  Evaluate compliance at each threshold of a density field.
%
%   result = sweepThresholdsReanalysis(x, analyses, thresholds, opts)
%
%   For each threshold t, constructs a binary density field (1 for solid
%   elements, xVoid for void), solves the FEM problem for every config,
%   and records volume fraction and compliance.  Optionally applies
%   removeDisconnectedComponents before each solve.
%
%   Inputs:
%     x           [nElems x 1]  full-arm density field (e.g. result.xFinal)
%     analyses    {nConfigs x 1} cell of FEAnalysis objects
%     thresholds  [nT x 1]  threshold values to evaluate (e.g. linspace(0.1,0.9,9))
%     opts        struct  (all fields optional)
%
%   opts fields:
%     pAgg         p-norm aggregation exponent              (default 4)
%     weights      [nConfigs x 1] per-config weights        (default uniform)
%     C0           [nConfigs x 1] compliance normalisation  (default: first threshold)
%     useParallel  logical, parallel config solves          (default false)
%     xVoid        void-element density for reanalysis      (default 1e-6)
%     mesh         struct with .elems — if provided, disconnected islands
%                  are removed before each reanalysis
%     anchorElems  [k x 1] required when mesh is provided
%
%   Output:  result struct with fields
%     thresholds   [nT x 1]  same as input
%     volFrac      [nT x 1]  fraction of solid elements after cleanup
%     C            [nT x nConfigs]  compliance for each threshold and config
%     J            [nT x 1]  normalised p-norm objective
%     C0           [nConfigs x 1]  normalisation used
%     bestIdx      scalar  index of threshold with lowest J
%
%   Note on computation time: each row of the sweep requires nConfigs FEM
%   solves.  For a 587k-DOF arm model with 6 configs, each row takes
%   ~1–5 min.  Use a coarse threshold grid (7–11 points) for interactive use.

    if nargin < 4, opts = struct(); end

    nConfigs = numel(analyses);
    pAgg        = optField(opts, 'pAgg',        4);
    weights     = optField(opts, 'weights',     ones(nConfigs, 1) / nConfigs);
    useParallel = optField(opts, 'useParallel', false);
    xVoid       = optField(opts, 'xVoid',       1e-6);
    doCleanup   = isfield(opts, 'mesh') && isfield(opts, 'anchorElems');

    weights = weights(:);
    thresholds = thresholds(:);
    nT = numel(thresholds);
    nElems = numel(x);

    volFrac = zeros(nT, 1);
    C       = zeros(nT, nConfigs);

    fprintf('sweepThresholdsReanalysis: %d thresholds × %d configs = %d FEM solves.\n', ...
        nT, nConfigs, nT * nConfigs);

    for i = 1:nT
        t     = thresholds(i);
        solid = x(:) >= t;

        if doCleanup
            solid = removeDisconnectedComponents(solid, opts.mesh, opts.anchorElems);
        end

        volFrac(i) = mean(solid);
        x_bin = double(solid) + (~solid) * xVoid;

        useParallelI = useParallel && license('test', 'Distrib_Computing_Toolbox');
        if useParallelI
            parfor k = 1:nConfigs
                C(i, k) = computeComplianceOnly(analyses{k}, x_bin, 1);
            end
        else
            for k = 1:nConfigs
                C(i, k) = computeComplianceOnly(analyses{k}, x_bin, 1);
            end
        end

        fprintf('  [%2d/%2d] t=%.3f  V=%.3f  C_max=%.4e\n', ...
            i, nT, t, volFrac(i), max(C(i, :)));
    end

    % Normalisation: use provided C0, or the compliance at the threshold
    % whose volume fraction is closest to the median volume fraction.
    C0 = optField(opts, 'C0', []);
    if isempty(C0)
        [~, refIdx] = min(abs(volFrac - median(volFrac)));
        C0 = max(C(refIdx, :)', eps);
        fprintf('  Normalising against threshold t=%.3f (median volume fraction).\n', ...
            thresholds(refIdx));
    end

    Cagg        = bsxfun(@rdivide, C, C0');       % [nT x nConfigs]
    weightedSum = (Cagg .^ pAgg) * weights;        % [nT x 1]
    J           = weightedSum .^ (1.0 / pAgg);

    [~, bestIdx] = min(J);

    result.thresholds = thresholds;
    result.volFrac    = volFrac;
    result.C          = C;
    result.J          = J;
    result.C0         = C0;
    result.bestIdx    = bestIdx;
    result.bestThreshold = thresholds(bestIdx);
    result.bestVolFrac   = volFrac(bestIdx);
end

% -------------------------------------------------------------------------
function v = optField(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

function result = solveSIMPComplianceVolumeMMA(analyses, z0, xmin, xmax, opts)
% solveSIMPComplianceVolumeMMA  Multi-config SIMP MMA: min compliance, vol.
%
%   z is the optimizer design vector.  For unlinked designs z is the full
%   density vector.  For linked designs opts.expand maps z -> full-arm x and
%   opts.pullback maps full-arm sensitivities back to z-space.

    opts = setDefaultOptions(opts, numel(z0));

    z = z0(:);
    xmin = xmin(:);
    xmax = xmax(:);
    n = numel(z);
    m = 1;
    fixedDesignVariables = normalizeFixedDesignVariables(opts.fixedDesignVariables, n);
    if ~isempty(fixedDesignVariables)
        xmin(fixedDesignVariables) = 1.0;
        xmax(fixedDesignVariables) = 1.0;
        z(fixedDesignVariables) = 1.0;
    end
    freeDesignVariables = setdiff((1:n)', fixedDesignVariables);
    nFree = numel(freeDesignVariables);
    assert(nFree > 0, 'No free design variables remain after applying fixedDesignVariables.');

    x0 = opts.expand(z);
    if isempty(opts.C0)
        [J0, dJdx0, C0, Cinit] = evaluateObjectiveAndGradient( ...
            analyses, x0, opts.penal, opts.pAgg, opts.weights, [], opts.useParallel);
    else
        [J0, dJdx0, C0, Cinit] = evaluateObjectiveAndGradient( ...
            analyses, x0, opts.penal, opts.pAgg, opts.weights, opts.C0, opts.useParallel);
    end
    grad0 = opts.pullback(dJdx0);

    zHistory = zeros(n, opts.maxIter + 1);
    zHistory(:, 1) = z;
    JHistory = nan(opts.maxIter + 1, 1);
    JHistory(1) = J0;
    CHistory = nan(opts.maxIter + 1, numel(analyses));
    CHistory(1, :) = Cinit(:)';
    volHistory = nan(opts.maxIter + 1, 1);
    volHistory(1) = mean(z);
    changeHistory = nan(opts.maxIter + 1, 1);
    changeHistory(1) = 0;
    iterTimeHistory = nan(opts.maxIter + 1, 1);
    iterTimeHistory(1) = 0;

    zFree = z(freeDesignVariables);
    xold1 = zFree;
    xold2 = zFree;
    low = zeros(nFree, 1);
    upp = ones(nFree, 1);
    a0 = 1;
    a = 0;
    c_mma = 1000;
    d = 0;
    objectiveScale = 1.0 / max(abs(J0), eps);

    nIter = 0;
    for iter = 1:opts.maxIter
        iterTic = tic;
        x = opts.expand(z);
        [J, dJdx] = evaluateObjectiveAndGradient( ...
            analyses, x, opts.penal, opts.pAgg, opts.weights, C0, opts.useParallel);
        gradJ = opts.sensitivityFilter * opts.pullback(dJdx);

        constr = sum(z) / (opts.VolFrac * n) - 1.0;
        gradConstr = ones(1, nFree) / (opts.VolFrac * n);

        zFree = z(freeDesignVariables);
        gradJFree = gradJ(freeDesignVariables);
        [zMmaFree, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = mmasub2( ...
            m, nFree, iter, zFree, xmin(freeDesignVariables), xmax(freeDesignVariables), xold1, xold2, ...
            objectiveScale * J, objectiveScale * gradJFree, 0 * gradJFree, ...
            constr, gradConstr, 0 * gradConstr, ...
            low, upp, a0, a, c_mma, d);

        if iter > 1
            xold2 = xold1;
        end
        xold1 = zFree;

        currentMoveLimit = max(opts.minMoveLimit, opts.moveLimit * opts.moveDecay^(iter - 1));
        zCandidateFree = min(max(zMmaFree, zFree - currentMoveLimit), zFree + currentMoveLimit);
        zCandidateFree = zFree + opts.mmaDamping * (zCandidateFree - zFree);
        zCandidate = z;
        zCandidate(freeDesignVariables) = zCandidateFree;
        zCandidate(fixedDesignVariables) = 1.0;
        zCandidate = enforceVolumeFraction(zCandidate, opts.VolFrac, xmin, xmax);
        zCandidate(fixedDesignVariables) = 1.0;

        change = max(abs(zCandidate - z));
        z = zCandidate;

        [Jnew, Cnew] = evaluateObjectiveOnly( ...
            analyses, opts.expand(z), opts.penal, opts.pAgg, opts.weights, C0, opts.useParallel);
        iterTimeSec = toc(iterTic);

        zHistory(:, iter + 1) = z;
        JHistory(iter + 1) = Jnew;
        CHistory(iter + 1, :) = Cnew(:)';
        volHistory(iter + 1) = mean(z);
        changeHistory(iter + 1) = change;
        iterTimeHistory(iter + 1) = iterTimeSec;

        printComplianceIteration(iter, Jnew, JHistory(iter), volHistory(iter + 1), ...
            change, iterTimeSec, Cnew, opts.configNames);

        nIter = iter;
        if iter >= opts.minIter && change < opts.changeTol
            fprintf('  Change tolerance %.1e reached at iteration %d (change=%.3e).\n', ...
                opts.changeTol, iter, change);
            break;
        end
    end

    finalIdx = nIter + 1;
    history.iteration = (0:nIter)';
    history.J = JHistory(1:finalIdx);
    history.C = CHistory(1:finalIdx, :);
    history.volumeFraction = volHistory(1:finalIdx);
    history.change = changeHistory(1:finalIdx);
    history.iterationTimeSec = iterTimeHistory(1:finalIdx);
    history.x = zHistory(:, 1:finalIdx);
    history.configNames = string(opts.configNames(:));
    history.configLabels = string(opts.configLabels(:));
    history.weights = opts.weights;

    result.zFinal = z;
    result.xFinal = opts.expand(z);
    result.finalJ = history.J(end);
    result.finalC = history.C(end, :)';
    result.finalVolumeFraction = mean(z);
    result.finalChange = history.change(end);
    result.nIter = nIter;
    result.history = history;
    result.C0 = C0;
    result.Cinit = Cinit;
    result.J0 = J0;
    result.grad0 = grad0;
    result.objectiveWindowConverged = objectiveHistoryConverged(history.J, opts.objectiveTol);
    result.objectiveDecreased = result.finalJ < history.J(1);
    result.mmaAcceptable = result.objectiveWindowConverged || result.objectiveDecreased;
end

function opts = setDefaultOptions(opts, n)
    if ~isfield(opts, 'expand') || isempty(opts.expand), opts.expand = @(z) z; end
    if ~isfield(opts, 'pullback') || isempty(opts.pullback), opts.pullback = @(dx) dx; end
    if ~isfield(opts, 'sensitivityFilter') || isempty(opts.sensitivityFilter)
        opts.sensitivityFilter = speye(n);
    end
    if ~isfield(opts, 'C0'), opts.C0 = []; end
    if ~isfield(opts, 'useParallel'), opts.useParallel = false; end
    if ~isfield(opts, 'minIter'), opts.minIter = 1; end
    if ~isfield(opts, 'objectiveTol'), opts.objectiveTol = 5.0e-3; end
    if ~isfield(opts, 'fixedDesignVariables'), opts.fixedDesignVariables = []; end
    if ~isfield(opts, 'configNames'), opts.configNames = strings(numel(opts.weights), 1); end
    if ~isfield(opts, 'configLabels'), opts.configLabels = opts.configNames; end
end

function fixed = normalizeFixedDesignVariables(fixed, n)
    fixed = unique(fixed(:));
    fixed = fixed(fixed >= 1 & fixed <= n);
end

function printComplianceIteration(iter, Jnew, Jold, vf, change, iterTimeSec, Cnew, configNames)
    fprintf('%4d  J=%12.6e  dJ/J0=% .3e  vf=%.4f  change=%.3e  time=%.2fs', ...
        iter, Jnew, (Jnew - Jold) / max(abs(Jold), eps), vf, change, iterTimeSec);
    for k = 1:numel(Cnew)
        fprintf('  C_%s=%.3e', char(configNames(k)), Cnew(k));
    end
    fprintf('\n');
end

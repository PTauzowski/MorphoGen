function result = solveSIMPVolumeComplianceMMA(analyses, z0, xmin, xmax, opts)
% solveSIMPVolumeComplianceMMA  Multi-config SIMP MMA: min volume, C limits.

    opts = setDefaultOptions(opts, numel(z0));

    z = z0(:);
    xmin = xmin(:);
    xmax = xmax(:);
    n = numel(z);
    m = numel(analyses);
    fixedDesignVariables = normalizeFixedDesignVariables(opts.fixedDesignVariables, n);
    if ~isempty(fixedDesignVariables)
        xmin(fixedDesignVariables) = 1.0;
        xmax(fixedDesignVariables) = 1.0;
        z(fixedDesignVariables) = 1.0;
    end
    freeDesignVariables = setdiff((1:n)', fixedDesignVariables);
    nFree = numel(freeDesignVariables);
    assert(nFree > 0, 'No free design variables remain after applying fixedDesignVariables.');

    [C0, dC0dx] = evaluateComplianceSetLocal(analyses, opts.expand(z), opts.penal, opts.useParallel);
    Ctarget = opts.C_coeff * max(C0, eps);
    grad0 = opts.pullback(dC0dx);
    objectiveScale = opts.objectiveScale;

    zHistory = zeros(n, opts.maxIter + 1);
    zHistory(:, 1) = z;
    volHistory = nan(opts.maxIter + 1, 1);
    volHistory(1) = mean(z);
    CHistory = nan(opts.maxIter + 1, m);
    CHistory(1, :) = C0(:)';
    constraintHistory = nan(opts.maxIter + 1, m);
    constraintHistory(1, :) = (C0(:) ./ Ctarget(:) - 1)';
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
    a = zeros(m, 1);
    c_mma = 1000 * ones(m, 1);
    d = zeros(m, 1);

    nIter = 0;
    for iter = 1:opts.maxIter
        iterTic = tic;
        [C, dCdx] = evaluateComplianceSetLocal(analyses, opts.expand(z), opts.penal, opts.useParallel);
        dCdz = opts.sensitivityFilter * opts.pullback(dCdx);

        objective = objectiveScale * mean(z);
        gradObjective = objectiveScale * ones(nFree, 1) / n;
        constr = C(:) ./ Ctarget(:) - 1;
        gradConstr = (dCdz(freeDesignVariables, :) ./ Ctarget(:)')';

        zFree = z(freeDesignVariables);
        [zMmaFree, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = mmasub2( ...
            m, nFree, iter, zFree, xmin(freeDesignVariables), xmax(freeDesignVariables), xold1, xold2, ...
            objective, gradObjective, 0 * gradObjective, ...
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
        zCandidate = min(max(zCandidate, xmin), xmax);
        zCandidate(fixedDesignVariables) = 1.0;

        change = max(abs(zCandidate - z));
        z = zCandidate;

        [Cnew, ~] = evaluateComplianceSetLocal(analyses, opts.expand(z), opts.penal, opts.useParallel);
        constrNew = Cnew(:) ./ Ctarget(:) - 1;
        iterTimeSec = toc(iterTic);

        zHistory(:, iter + 1) = z;
        volHistory(iter + 1) = mean(z);
        CHistory(iter + 1, :) = Cnew(:)';
        constraintHistory(iter + 1, :) = constrNew(:)';
        changeHistory(iter + 1) = change;
        iterTimeHistory(iter + 1) = iterTimeSec;
        nIter = iter;

        printVolumeIteration(iter, volHistory(iter + 1), constrNew, change, ...
            iterTimeSec, Cnew, opts.configNames);

        activeEnough = max(constrNew) >= -opts.activeConstraintTol;
        if iter >= opts.minIter && change < opts.changeTol && ...
                max(constrNew) <= opts.constraintTol && activeEnough
            fprintf(['  Change tolerance %.1e reached with feasible and active constraints ' ...
                'at iteration %d.\n'], opts.changeTol, iter);
            break;
        end
    end

    finalIdx = nIter + 1;
    history.iteration = (0:nIter)';
    history.volumeFraction = volHistory(1:finalIdx);
    history.C = CHistory(1:finalIdx, :);
    history.C0 = C0;
    history.Ctarget = Ctarget;
    history.constraint = constraintHistory(1:finalIdx, :);
    history.change = changeHistory(1:finalIdx);
    history.iterationTimeSec = iterTimeHistory(1:finalIdx);
    history.x = zHistory(:, 1:finalIdx);
    history.configNames = string(opts.configNames(:));
    history.configLabels = string(opts.configLabels(:));

    result.zFinal = z;
    result.xFinal = opts.expand(z);
    result.finalVolumeFraction = mean(z);
    result.finalC = history.C(end, :)';
    result.finalConstraint = history.constraint(end, :)';
    result.finalChange = history.change(end);
    result.constraintsSatisfied = max(result.finalConstraint) <= opts.constraintTol;
    result.nIter = nIter;
    result.history = history;
    result.C0 = C0;
    result.Ctarget = Ctarget;
    result.grad0 = grad0;
end

function opts = setDefaultOptions(opts, n)
    if ~isfield(opts, 'expand') || isempty(opts.expand), opts.expand = @(z) z; end
    if ~isfield(opts, 'pullback') || isempty(opts.pullback), opts.pullback = @(dx) dx; end
    if ~isfield(opts, 'sensitivityFilter') || isempty(opts.sensitivityFilter)
        opts.sensitivityFilter = speye(n);
    end
    if ~isfield(opts, 'useParallel'), opts.useParallel = false; end
    if ~isfield(opts, 'minIter'), opts.minIter = 1; end
    if ~isfield(opts, 'constraintTol'), opts.constraintTol = 5.0e-3; end
    if ~isfield(opts, 'activeConstraintTol'), opts.activeConstraintTol = 2.0e-2; end
    if ~isfield(opts, 'objectiveScale'), opts.objectiveScale = n; end
    if ~isfield(opts, 'fixedDesignVariables'), opts.fixedDesignVariables = []; end
    if ~isfield(opts, 'configNames'), opts.configNames = strings(0, 1); end
    if ~isfield(opts, 'configLabels'), opts.configLabels = opts.configNames; end
end

function fixed = normalizeFixedDesignVariables(fixed, n)
    fixed = unique(fixed(:));
    fixed = fixed(fixed >= 1 & fixed <= n);
end

function [C, dC] = evaluateComplianceSetLocal(analyses, x, penal, useParallel)
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
end

function printVolumeIteration(iter, vf, constrNew, change, iterTimeSec, Cnew, configNames)
    fprintf('%4d  vf=%.4f  maxConstr=% .3e  change=%.3e  time=%.2fs', ...
        iter, vf, max(constrNew), change, iterTimeSec);
    for k = 1:numel(Cnew)
        fprintf('  C_%s=%.3e', char(configNames(k)), Cnew(k));
    end
    fprintf('\n');
end

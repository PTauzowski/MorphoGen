function result = solveSIMPVolumeStressMMA(analyses, z0, xmin, xmax, opts)
% solveSIMPVolumeStressMMA  Multi-config SIMP MMA: min volume, stress limits.
%
%   Minimize volume fraction subject to P-norm stress aggregates <= target.
%   Stress sensitivities are computed via the adjoint method (exact).

    opts = setDefaultOptions(opts, numel(z0));

    z = z0(:);
    xmin = xmin(:);
    xmax = xmax(:);
    n = numel(z);
    nConfigs = numel(analyses);
    nStressClusters = opts.nStressClusters;

    fixedDesignVariables = normalizeFixedDesignVariables(opts.fixedDesignVariables, n);
    if ~isempty(fixedDesignVariables)
        xmin(fixedDesignVariables) = 1.0;
        xmax(fixedDesignVariables) = 1.0;
        z(fixedDesignVariables) = 1.0;
    end
    freeDesignVariables = setdiff((1:n)', fixedDesignVariables);
    nFree = numel(freeDesignVariables);
    assert(nFree > 0, 'No free design variables remain after applying fixedDesignVariables.');

    [zPhysical, dzPhysicalDz] = physicalDesign(z, opts, fixedDesignVariables);
    x0 = opts.expand(zPhysical);
    [S0, dS0dx, maxStress0] = evaluateStressSetLocal( ...
        analyses, x0, opts.penal, opts.stressPNorm, opts.stressRelaxationQ, ...
        [], opts.useParallel, nStressClusters);
    m = numel(S0);
    opts.constraintNames = stressConstraintNames(opts.configNames, nStressClusters);
    Starget = stressTargets(S0, opts.stressCoeff, nConfigs, nStressClusters);
    grad0 = chainDesignGradient(dS0dx, opts, dzPhysicalDz);

    zHistory = zeros(n, opts.maxIter + 1);
    zHistory(:, 1) = z;
    zPhysicalHistory = zeros(n, opts.maxIter + 1);
    zPhysicalHistory(:, 1) = zPhysical;
    volHistory = nan(opts.maxIter + 1, 1);
    volHistory(1) = mean(zPhysical);
    SHistory = nan(opts.maxIter + 1, m);
    SHistory(1, :) = S0(:)';
    maxStressHistory = nan(opts.maxIter + 1, m);
    maxStressHistory(1, :) = maxStress0(:)';
    constraintHistory = nan(opts.maxIter + 1, m);
    constraintHistory(1, :) = (S0(:) ./ Starget(:) - 1)';
    changeHistory = nan(opts.maxIter + 1, 1);
    changeHistory(1) = 0;
    iterTimeHistory = nan(opts.maxIter + 1, 1);
    iterTimeHistory(1) = 0;

    % Best feasible design tracking
    bestFeasible.z = z;
    bestFeasible.zPhysical = zPhysical;
    bestFeasible.volFrac = mean(zPhysical);
    bestFeasible.iter = 0;
    bestFeasible.valid = max(S0(:) ./ Starget(:) - 1) <= opts.constraintTol;

    zFree = z(freeDesignVariables);
    xold1 = zFree;
    xold2 = zFree;
    low = zeros(nFree, 1);
    upp = ones(nFree, 1);
    a0 = 1;
    a = zeros(m, 1);
    c_mma = opts.mmaConstraintScale * ones(m, 1);
    d = zeros(m, 1);
    objectiveScale = opts.objectiveScale;

    nIter = 0;
    for iter = 1:opts.maxIter
        iterTic = tic;

        [zPhysical, dzPhysicalDz] = physicalDesign(z, opts, fixedDesignVariables);
        x = opts.expand(zPhysical);
        [S, dSdx, ~] = evaluateStressSetLocal( ...
            analyses, x, opts.penal, opts.stressPNorm, opts.stressRelaxationQ, ...
            Starget, opts.useParallel, nStressClusters);
        dSdz = chainDesignGradient(dSdx, opts, dzPhysicalDz);

        objective = objectiveScale * mean(zPhysical);
        gradObjectiveAll = chainDesignSpaceGradient(ones(n, 1) / n, opts, dzPhysicalDz);
        gradObjective = objectiveScale * gradObjectiveAll(freeDesignVariables);
        constr = S(:) ./ Starget(:) - 1;
        gradConstr = (dSdz(freeDesignVariables, :) ./ Starget(:)')';

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

        [zPhysicalNew, ~] = physicalDesign(z, opts, fixedDesignVariables);
        [Snew, ~, maxStressNew] = evaluateStressSetLocal( ...
            analyses, opts.expand(zPhysicalNew), opts.penal, opts.stressPNorm, ...
            opts.stressRelaxationQ, Starget, opts.useParallel, nStressClusters);
        constrNew = Snew(:) ./ Starget(:) - 1;
        iterTimeSec = toc(iterTic);

        zHistory(:, iter + 1) = z;
        zPhysicalHistory(:, iter + 1) = zPhysicalNew;
        volHistory(iter + 1) = mean(zPhysicalNew);
        SHistory(iter + 1, :) = Snew(:)';
        maxStressHistory(iter + 1, :) = maxStressNew(:)';
        constraintHistory(iter + 1, :) = constrNew(:)';
        changeHistory(iter + 1) = change;
        iterTimeHistory(iter + 1) = iterTimeSec;
        nIter = iter;

        printStressIteration(iter, volHistory(iter + 1), constrNew, change, ...
            iterTimeSec, Snew, maxStressNew, opts.constraintNames);

        % Track best feasible design
        if max(constrNew) <= opts.constraintTol
            if ~bestFeasible.valid || mean(zPhysicalNew) < bestFeasible.volFrac
                bestFeasible.z = z;
                bestFeasible.zPhysical = zPhysicalNew;
                bestFeasible.volFrac = mean(zPhysicalNew);
                bestFeasible.iter = iter;
                bestFeasible.valid = true;
            end
        end

        % Convergence check
        activeEnough = max(constrNew) >= -opts.activeConstraintTol;
        if iter >= opts.minIter && change < opts.changeTol && ...
                max(constrNew) <= opts.constraintTol && activeEnough
            fprintf(['  Change tolerance %.1e reached with feasible and active stress ' ...
                'constraints at iteration %d.\n'], opts.changeTol, iter);
            break;
        end

        % Divergence guard: stop if constraint max has been rising for divergenceWindow iters
        if iter >= opts.minIter + opts.divergenceWindow
            windowStart = iter - opts.divergenceWindow + 2;
            windowEnd   = iter + 1;
            recentMaxConstr = max(constraintHistory(windowStart:windowEnd, :), [], 2);
            if recentMaxConstr(end) > opts.divergenceTol && all(diff(recentMaxConstr) > 0)
                fprintf(['  Divergence detected: max stress constraint has increased for ' ...
                    '%d consecutive iterations (last=%.3e). Stopping.\n'], ...
                    opts.divergenceWindow, recentMaxConstr(end));
                break;
            end
        end
    end

    finalIdx = nIter + 1;
    history.iteration = (0:nIter)';
    history.volumeFraction = volHistory(1:finalIdx);
    history.stressAggregate = SHistory(1:finalIdx, :);
    history.initialStressAggregate = S0;
    history.targetStressAggregate = Starget;
    history.maxStress = maxStressHistory(1:finalIdx, :);
    history.constraint = constraintHistory(1:finalIdx, :);
    history.change = changeHistory(1:finalIdx);
    history.iterationTimeSec = iterTimeHistory(1:finalIdx);
    history.x = zHistory(:, 1:finalIdx);
    history.xPhysical = zPhysicalHistory(:, 1:finalIdx);
    history.configNames = string(opts.constraintNames(:));
    history.configLabels = string(opts.constraintNames(:));
    history.loadConfigNames = string(opts.configNames(:));
    history.loadConfigLabels = string(opts.configLabels(:));
    history.nStressClusters = nStressClusters;

    [zPhysicalFinal, ~] = physicalDesign(z, opts, fixedDesignVariables);
    result.zFinal = z;
    result.zPhysicalFinal = zPhysicalFinal;
    result.xFinal = opts.expand(zPhysicalFinal);
    result.finalVolumeFraction = mean(zPhysicalFinal);
    result.finalStressAggregate = history.stressAggregate(end, :)';
    result.finalMaxStress = history.maxStress(end, :)';
    result.finalConstraint = history.constraint(end, :)';
    result.finalChange = history.change(end);
    result.constraintsSatisfied = max(result.finalConstraint) <= opts.constraintTol;
    result.nIter = nIter;
    result.history = history;
    result.S0 = S0;
    result.Starget = Starget;
    result.maxStress0 = maxStress0;
    result.grad0 = grad0;
    result.bestFeasible = bestFeasible;
    result.projection = struct( ...
        'enabled', opts.useProjection, ...
        'beta', opts.projectionBeta, ...
        'eta', opts.projectionEta);
    result.nStressClusters = nStressClusters;
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
    if ~isfield(opts, 'stressPNorm'), opts.stressPNorm = 12; end
    if ~isfield(opts, 'stressRelaxationQ'), opts.stressRelaxationQ = 0.5; end
    if ~isfield(opts, 'stressCoeff'), opts.stressCoeff = 2.0; end
    if ~isfield(opts, 'mmaConstraintScale'), opts.mmaConstraintScale = 1000; end
    if ~isfield(opts, 'configNames'), opts.configNames = strings(0, 1); end
    if ~isfield(opts, 'configLabels'), opts.configLabels = opts.configNames; end
    if ~isfield(opts, 'useProjection'), opts.useProjection = false; end
    if ~isfield(opts, 'projectionBeta'), opts.projectionBeta = 1.0; end
    if ~isfield(opts, 'projectionEta'), opts.projectionEta = 0.5; end
    if ~isfield(opts, 'nStressClusters'), opts.nStressClusters = 1; end
    if ~isfield(opts, 'divergenceTol'), opts.divergenceTol = 0.15; end
    if ~isfield(opts, 'divergenceWindow'), opts.divergenceWindow = 20; end
    opts.nStressClusters = max(1, round(opts.nStressClusters));
end

function fixed = normalizeFixedDesignVariables(fixed, n)
    fixed = unique(fixed(:));
    fixed = fixed(fixed >= 1 & fixed <= n);
end

function [zPhysical, dzPhysicalDz] = physicalDesign(z, opts, fixedDesignVariables)
    z = z(:);
    if opts.useProjection
        zFiltered = opts.sensitivityFilter * z;
        zFiltered = min(max(zFiltered, 0.0), 1.0);
        [zPhysical, dzPhysicalDz] = heavisideProjection( ...
            zFiltered, opts.projectionBeta, opts.projectionEta);
        if ~isempty(fixedDesignVariables)
            zPhysical(fixedDesignVariables) = 1.0;
            dzPhysicalDz(fixedDesignVariables) = 0.0;
        end
    else
        zPhysical = z;
        dzPhysicalDz = ones(size(z));
    end
end

function [rho, drho] = heavisideProjection(rhoTilde, beta, eta)
    beta = max(beta, eps);
    eta = min(max(eta, eps), 1.0 - eps);
    denom = tanh(beta * eta) + tanh(beta * (1.0 - eta));
    rho = (tanh(beta * eta) + tanh(beta * (rhoTilde - eta))) / denom;
    drho = beta * (1.0 - tanh(beta * (rhoTilde - eta)).^2) / denom;
end

function gradZ = chainDesignGradient(gradPhysicalFull, opts, dzPhysicalDz)
    gradPhysical = opts.pullback(gradPhysicalFull);
    if opts.useProjection
        gradZ = opts.sensitivityFilter' * (dzPhysicalDz .* gradPhysical);
    else
        gradZ = opts.sensitivityFilter * gradPhysical;
    end
end

function gradZ = chainDesignSpaceGradient(gradPhysicalDesign, opts, dzPhysicalDz)
    if opts.useProjection
        gradZ = opts.sensitivityFilter' * (dzPhysicalDz .* gradPhysicalDesign);
    else
        gradZ = gradPhysicalDesign;
    end
end

function names = stressConstraintNames(configNames, nStressClusters)
    configNames = string(configNames(:));
    if nStressClusters == 1
        names = configNames;
        return;
    end

    names = strings(numel(configNames) * nStressClusters, 1);
    c = 0;
    for k = 1:numel(configNames)
        for j = 1:nStressClusters
            c = c + 1;
            names(c) = sprintf('%s_c%02d', char(configNames(k)), j);
        end
    end
end

function Starget = stressTargets(S0, stressCoeff, nConfigs, nStressClusters)
    S0 = S0(:);
    Starget = zeros(size(S0));
    for k = 1:nConfigs
        idx = (k - 1) * nStressClusters + (1:nStressClusters);
        target = stressCoeff * max(S0(idx));
        Starget(idx) = max(target, eps);
    end
end

function [S, dS, maxStress] = evaluateStressSetLocal( ...
        analyses, x, penal, stressPNorm, stressRelaxationQ, Starget, useParallel, nStressClusters)
    nConfigs = numel(analyses);
    nConstraints = nConfigs * nStressClusters;
    S = zeros(nConstraints, 1);
    dS = zeros(numel(x), nConstraints);
    maxStress = zeros(nConstraints, 1);

    useParallel = useParallel && license('test', 'Distrib_Computing_Toolbox');
    if useParallel
        SCell = cell(nConfigs, 1);
        dSCell = cell(nConfigs, 1);
        maxStressCell = cell(nConfigs, 1);
        parfor k = 1:nConfigs
            idx = (k - 1) * nStressClusters + (1:nStressClusters);
            targetK = [];
            if ~isempty(Starget)
                targetK = Starget(idx);
            end
            [SCell{k}, dSCell{k}, maxStressCell{k}] = computeStressAggregateAndGradient( ...
                analyses{k}, x, penal, stressPNorm, stressRelaxationQ, targetK, nStressClusters);
        end
        for k = 1:nConfigs
            idx = (k - 1) * nStressClusters + (1:nStressClusters);
            S(idx) = SCell{k};
            dS(:, idx) = dSCell{k};
            maxStress(idx) = maxStressCell{k};
        end
    else
        for k = 1:nConfigs
            idx = (k - 1) * nStressClusters + (1:nStressClusters);
            targetK = [];
            if ~isempty(Starget)
                targetK = Starget(idx);
            end
            [S(idx), dS(:, idx), maxStress(idx)] = computeStressAggregateAndGradient( ...
                analyses{k}, x, penal, stressPNorm, stressRelaxationQ, targetK, nStressClusters);
        end
    end
end

function [S, dS, maxStress] = computeStressAggregateAndGradient( ...
        analysis, x, penal, stressPNorm, q, Starget, nStressClusters)
    nElems = analysis.getTotalElemsNumber();
    assert(numel(x) == nElems, ...
        'Density vector length %d does not match analysis element count %d.', ...
        numel(x), nElems);

    x = x(:);
    xPenal = x .^ penal;
    analysis.solveWeighted(xPenal, true);
    analysis.computeElementResults(xPenal);

    % Mean-GP Huber-Mises stress of stored (penalized) stress
    sigma = elementGPHuberMisesStress(analysis, nElems);

    % qp-relaxed stress: x^q * sigma_stored_HM = x^{p+q} * sigma_phys_HM
    relaxed = (max(x, eps) .^ q) .* sigma;

    if isempty(Starget)
        target = ones(nStressClusters, 1);
    else
        target = max(Starget(:), eps);
    end

    clusters = stressLevelClusters(relaxed, nStressClusters);
    S = zeros(nStressClusters, 1);
    dSdRelaxedAll = zeros(nElems, nStressClusters);
    maxStress = zeros(nStressClusters, 1);

    for c = 1:nStressClusters
        elemIds = clusters{c};
        maxStress(c) = max(sigma(elemIds));
        ratio = relaxed(elemIds) / target(c);
        meanPower = mean(ratio .^ stressPNorm);
        S(c) = target(c) * meanPower ^ (1.0 / stressPNorm);

        if meanPower <= eps
            continue;
        end

        dSdr = target(c) * (1.0 / stressPNorm) * meanPower^(1.0 / stressPNorm - 1.0) ...
            * stressPNorm * ratio.^(stressPNorm - 1.0) / numel(elemIds);
        dSdRelaxedAll(elemIds, c) = dSdr / target(c);
    end

    % Direct term: d(relaxed_e)/dx_e = (p+q)/x_e * relaxed_e (holding u fixed)
    % relaxed = x^q * sigma_stored = x^{p+q} * sigma_phys, so d/dx = (p+q)/x * relaxed
    dS_direct = dSdRelaxedAll .* ((penal + q) * relaxed ./ max(x, eps));

    % Indirect term via adjoint: d(sigma)/du contributes through K*u=P
    P_adj = assembleStressAdjointLoad(analysis, x, penal, q, dSdRelaxedAll, sigma);
    lambda_fem = analysis.solveAdjointWithLoad(xPenal, P_adj);
    dS_indirect = computeAdjointElementSensitivity(analysis, lambda_fem, x, penal);

    dS = dS_direct + dS_indirect;
end

function clusters = stressLevelClusters(stressMeasure, nStressClusters)
    [~, order] = sort(stressMeasure(:), 'descend');
    nElems = numel(order);
    nStressClusters = min(max(1, nStressClusters), nElems);
    clusters = cell(nStressClusters, 1);
    edges = round(linspace(0, nElems, nStressClusters + 1));
    for c = 1:nStressClusters
        ids = order((edges(c) + 1):edges(c + 1));
        if isempty(ids)
            ids = order(1);
        end
        clusters{c} = ids;
    end
end

function sigma = elementGPHuberMisesStress(analysis, nElems)
    sigma = zeros(nElems, 1);
    elemIndices = analysis.getElemIndices();
    for i = 1:numel(analysis.felems)
        fe = analysis.felems{i};
        if ~isfield(fe.results, 'gp') || ~isfield(fe.results.gp, 'stress')
            continue;
        end
        s = fe.results.gp.stress;   % (nElemsI, nip, 6)
        s1 = s(:,:,1); s2 = s(:,:,2); s3 = s(:,:,3);
        s4 = s(:,:,4); s5 = s(:,:,5); s6 = s(:,:,6);
        hmGP = sqrt(0.5*((s1-s2).^2 + (s2-s3).^2 + (s3-s1).^2) + 3*(s4.^2+s5.^2+s6.^2));
        elemIds = elemIndices{i};
        sigma(elemIds) = mean(hmGP, 2);   % mean over integration points
    end
end

function P_adj = assembleStressAdjointLoad(analysis, x, penal, q, dSdRelaxedAll, ~)
    % Exact adjoint load for sigma_e = (1/nip)*sum_ip sHM(s_ip).
    % Uses per-IP v_ip = grad_sigma(sHM(s_ip)) so that
    % P_adj_c(j) = sum_e wc(e) * (1/nip)*sum_ip B_ip^T * D * v_ip(e).
    nTaskDim  = analysis.getTaskDim();
    nClusters = size(dSdRelaxedAll, 2);
    ndofsNode = size(analysis.ndofs, 2);
    P_adj = zeros(nTaskDim, nClusters);

    elemIndices = analysis.getElemIndices();
    nodes = analysis.mesh.nodes;

    for fi = 1:numel(analysis.felems)
        fe = analysis.felems{fi};
        if ~isa(fe, 'SolidElasticElem')
            continue;
        end

        elemIds = elemIndices{fi};
        nElemsI = numel(elemIds);
        nnpE    = size(fe.elems, 2);
        dimE    = ndofsNode * nnpE;

        s = fe.results.gp.stress;   % (nElemsI, nip, 6)
        D = fe.mat.D;

        % B-matrix quadrature setup
        integrator = fe.sf.createIntegrator();
        dN   = fe.sf.computeGradient(integrator.points);
        dNtr = permute(dN, [2, 1, 3]);
        nip  = size(integrator.points, 1);

        % f_all(dimE, nElemsI) = (1/nip)*sum_ip B_ip^T * D * v_ip per element.
        f_all   = zeros(dimE, nElemsI);
        chunkSz = fe.assemblyChunkSize(dimE);

        for first = 1:chunkSz:nElemsI
            chunk = first:min(first + chunkSz - 1, nElemsI);
            nc = numel(chunk);
            [Jinv, ~] = fe.jacobianInversePages(nodes, chunk, dNtr);
            B = fe.strainBPages(Jinv, dNtr);  % (6, dimE, nc, nip)

            for ip = 1:nip
                % Per-IP HM gradient direction for each element in chunk
                sip = reshape(s(chunk, ip, :), nc, 6);
                s1 = sip(:,1); s2 = sip(:,2); s3 = sip(:,3);
                s4 = sip(:,4); s5 = sip(:,5); s6 = sip(:,6);
                sHMip = sqrt(0.5*((s1-s2).^2+(s2-s3).^2+(s3-s1).^2)+3*(s4.^2+s5.^2+s6.^2));
                vip = [(2*s1-s2-s3),(2*s2-s1-s3),(2*s3-s1-s2),6*s4,6*s5,6*s6] ...
                      ./ (2 * max(sHMip, eps));              % (nc, 6)
                DvipPages = reshape(D * vip', 6, 1, nc);    % (6, 1, nc)

                Bip  = B(:, :, :, ip);                      % (6, dimE, nc)
                BtDv = reshape(squeeze(pagemtimes(Bip, 'transpose', DvipPages, 'none')), dimE, nc);
                f_all(:, chunk) = f_all(:, chunk) + BtDv / nip;
            end
        end

        % Global DOF map: allGlobalDofs(e, ld) = FEM DOF index
        nodeIdxPerDof = ceil((1:dimE)' / ndofsNode);   % (dimE,)
        dofIdxPerDof  = mod((0:dimE-1)', ndofsNode) + 1;
        allGlobalDofs = (fe.elems(:, nodeIdxPerDof') - 1) * ndofsNode + dofIdxPerDof';
        % allGlobalDofs: (nElemsI, dimE)

        % Weight and accumulate for each cluster
        xpq = x(elemIds) .^ (penal + q);   % (nElemsI, 1)
        for c = 1:nClusters
            wc = dSdRelaxedAll(elemIds, c) .* xpq;   % (nElemsI, 1)
            if all(wc == 0), continue; end
            weighted_f = f_all .* wc';                % (dimE, nElemsI)
            % allGlobalDofs is (nElemsI x dimE), weighted_f is (dimE x nElemsI).
            % Column-major (:) on allGlobalDofs visits (ld varies slowest, e fastest).
            % Transpose weighted_f before flattening to match that same ordering.
            P_adj(:, c) = P_adj(:, c) + accumarray(allGlobalDofs(:), reshape(weighted_f', [], 1), [nTaskDim, 1]);
        end
    end
end

function dS_indirect = computeAdjointElementSensitivity(analysis, lambda_fem, x, penal)
    % Indirect term: dS_c/dx_e = -penal * x_e^{p-1} * lambda_c^T * K_e^0 * u_e
    nElems   = numel(x);
    nClusters = size(lambda_fem, 2);
    dS_indirect = zeros(nElems, nClusters);

    elemIndices = analysis.getElemIndices();
    nodes   = analysis.mesh.nodes;
    u_nodal = analysis.qnodal;   % (nnodes, ndofs) — forward solution

    stiffnessFunction = 'computeStifnessMatrix';
    if analysis.isConst
        stiffnessFunction = 'computeStifnessMatrixConst';
    end

    for fi = 1:numel(analysis.felems)
        fe = analysis.felems{fi};
        if ~isa(fe, 'SolidElasticElem')
            continue;
        end

        elemIds = elemIndices{fi};
        nElemsI = numel(elemIds);
        nnpE    = size(fe.elems, 2);
        ndofsNode = numel(fe.ndofs);
        dim     = nnpE * ndofsNode;

        % Un-penalized element stiffness (x=1)
        K0flat = fe.(stiffnessFunction)(nodes, ones(nElemsI, 1));
        K0 = reshape(K0flat, dim, dim, nElemsI);

        % Forward displacement for this FE group
        uElems = fe.createElemSolutionVectors(u_nodal);   % (dim, nElemsI)
        K0u = pagemtimes(K0, reshape(uElems, dim, 1, nElemsI));  % (dim, 1, nElemsI)

        prefactor = -penal * x(elemIds) .^ (penal - 1);   % (nElemsI,)

        for c = 1:nClusters
            lambda_nodal_c = analysis.fromFEMVector(lambda_fem(:, c));  % (nnodes, ndofs)
            lambdaElems = fe.createElemSolutionVectors(lambda_nodal_c); % (dim, nElemsI)
            bilinear = squeeze(pagemtimes(reshape(lambdaElems, 1, dim, nElemsI), K0u));
            dS_indirect(elemIds, c) = dS_indirect(elemIds, c) + prefactor .* bilinear(:);
        end
    end
end

function printStressIteration(iter, vf, constrNew, change, iterTimeSec, Snew, maxStressNew, configNames)
    fprintf('%4d  vf=%.4f  maxStressConstr=% .3e  change=%.3e  time=%.2fs', ...
        iter, vf, max(constrNew), change, iterTimeSec);
    if numel(Snew) > 12
        [~, activeIdx] = max(constrNew);
        fprintf('  active=%s  S_active=%.3e  HMmax_active=%.3e', ...
            char(configNames(activeIdx)), Snew(activeIdx), maxStressNew(activeIdx));
        fprintf('  HMmax_global=%.3e', max(maxStressNew));
        fprintf('\n');
        return;
    end
    for k = 1:numel(Snew)
        fprintf('  S_%s=%.3e  HMmax_%s=%.3e', ...
            char(configNames(k)), Snew(k), char(configNames(k)), maxStressNew(k));
    end
    fprintf('\n');
end

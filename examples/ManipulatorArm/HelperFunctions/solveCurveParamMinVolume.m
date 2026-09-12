function [pBest, vfBest, historyRows] = solveCurveParamMinVolume( ...
    analyses, modelRef, curveOpts, stressLimit, dispLimit, ...
    paramBounds, p0, solverOpts)
% solveCurveParamMinVolume  Minimize volume fraction via SQP subject to
%   per-configuration stress and displacement constraints.
%
%   [pBest, vfBest, historyRows] = solveCurveParamMinVolume(
%       analyses, modelRef, curveOpts, stressLimit, dispLimit,
%       paramBounds, p0, solverOpts)
%
%   Inputs
%     analyses    - cell array of LinearElasticityWeighted, one per config
%     modelRef    - ManipulatorModel3D reference (first config)
%     curveOpts   - struct passed to buildCurveLinkedDensity /
%                   evaluateLinkedDensityMetrics; VolFrac field is updated
%                   from the optimized theta(end) at each evaluation
%     stressLimit - scalar stress upper bound [Pa]
%     dispLimit   - scalar tip-displacement upper bound [m] (Inf = inactive)
%     paramBounds - struct from defaultCurveParamBounds
%     p0          - initial curve parameter struct
%     solverOpts  - struct with fields:
%                     maxIter      (default 80)
%                     maxFunEvals  (default 400)
%                     vfInitial    (default 1.0)
%                     vfMin        (default 0.05)
%
%   Outputs
%     pBest       - optimal curve parameter struct
%     vfBest      - optimal volume fraction
%     historyRows - struct array with one row per function evaluation

    if nargin < 8 || isempty(solverOpts)
        solverOpts = struct();
    end
    maxIter     = localOpt(solverOpts, 'maxIter',     80);
    maxFunEvals = localOpt(solverOpts, 'maxFunEvals', 400);
    vfInitial   = localOpt(solverOpts, 'vfInitial',   1.0);
    vfMin       = localOpt(solverOpts, 'vfMin',       0.05);

    u0 = curveParamsToUnconstrained(p0, paramBounds);
    theta0 = [u0(:); vfInitial];
    lb = [-Inf(numel(u0), 1); vfMin];
    ub = [ Inf(numel(u0), 1); 1.0  ];

    % ratios in log are relative to the constraint limits (>1 = violated)

    evalCount   = 0;
    historyRows = struct([]);
    cachedTheta   = [];
    cachedMetrics = [];
    lastFeasibleTheta  = [];   % best feasible iterate seen so far
    lastFeasibleVf     = Inf;

    fminconOpts = optimoptions('fmincon', ...
        'Algorithm',                'sqp', ...
        'Display',                  'iter', ...
        'MaxIterations',            maxIter, ...
        'MaxFunctionEvaluations',   maxFunEvals, ...
        'OptimalityTolerance',      1e-2, ...
        'ConstraintTolerance',      3e-2, ...
        'FiniteDifferenceStepSize', 0.02, ...
        'FiniteDifferenceType',     'forward');

    fprintf('\nRunning fmincon (SQP): minimize VF, stressLimit=%.4e Pa\n', stressLimit);

    if maxIter <= 0 || maxFunEvals <= 1
        fprintf('Skipping fmincon (maxIter=%d, maxFunEvals=%d).\n', maxIter, maxFunEvals);
        thetaBest = theta0;
    else
        [thetaSqp, ~, exitflag] = fmincon(@objFn, theta0, [], [], [], [], lb, ub, ...
            @conFn, fminconOpts);
        % exitflag < 0  →  infeasible termination; fall back to last feasible.
        if exitflag < 0 && ~isempty(lastFeasibleTheta)
            fprintf('  fmincon exited infeasibly (flag %d); using last feasible iterate.\n', exitflag);
            thetaBest = lastFeasibleTheta;
        else
            thetaBest = thetaSqp;
        end
    end

    uBest  = thetaBest(1:end-1);
    vfBest = thetaBest(end);
    pBest  = unconstrainedToCurveParams(uBest, paramBounds);

    % ---- nested: objective ------------------------------------------------
    function J = objFn(theta)
        J = theta(end);
    end

    % ---- nested: constraints ----------------------------------------------
    function [c, ceq] = conFn(theta)
        m = evalCached(theta);
        c = m.maxHM(:) / stressLimit - 1;
        if isfinite(dispLimit)
            c = [c; abs(m.tipUz(:)) / dispLimit - 1];
        end
        ceq = [];
    end

    % ---- nested: cached FEM evaluation ------------------------------------
    function m = evalCached(theta)
        if ~isequal(theta, cachedTheta)
            params = unconstrainedToCurveParams(theta(1:end-1), paramBounds);
            vf     = max(vfMin, min(1.0, theta(end)));
            localCurveOpts         = curveOpts;
            localCurveOpts.VolFrac = vf;
            [~, rhoFull, info] = buildCurveLinkedDensity(params, modelRef, localCurveOpts);
            [metrics, ~]       = evaluateLinkedDensityMetrics(analyses, rhoFull, localCurveOpts);
            cachedTheta   = theta;
            cachedMetrics = metrics;

            evalCount = evalCount + 1;
            sRatio = max(metrics.maxHM)      / stressLimit;
            dRatio = max(abs(metrics.tipUz)) / max(abs(dispLimit), eps);

            % Track last feasible iterate for fallback on infeasible exit.
            vfNow = max(vfMin, min(1.0, theta(end)));
            isFeasible = all(metrics.maxHM(:) <= stressLimit) && ...
                (~isfinite(dispLimit) || all(abs(metrics.tipUz(:)) <= dispLimit));
            if isFeasible && vfNow < lastFeasibleVf
                lastFeasibleTheta = theta;
                lastFeasibleVf    = vfNow;
            end
            fprintf('  eval %4d: VF=%.4f  s_ratio=%.3f  d_ratio=%.3f\n', ...
                evalCount, info.fullVolumeFraction, sRatio, dRatio);

            row             = params;
            row.eval        = evalCount;
            row.VF          = info.fullVolumeFraction;
            row.stressRatio = sRatio;
            row.dispRatio   = dRatio;
            row.maxHM       = max(metrics.maxHM);
            row.maxAbsUz    = max(abs(metrics.tipUz));
            historyRows = [historyRows; row]; %#ok<AGROW>
        end
        m = cachedMetrics;
    end

end

function v = localOpt(s, name, default)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

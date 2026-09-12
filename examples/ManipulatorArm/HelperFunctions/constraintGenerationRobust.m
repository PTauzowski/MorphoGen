function [pBest, vfBest, adversarialConfigs, historyRows] = constraintGenerationRobust( ...
    modelRef, arm, curveOpts, stressLimit, dispLimit, ...
    paramBounds, p0, deltaDeg, cgOpts)
% constraintGenerationRobust  Robust topology optimization via constraint generation.
%
%   [pBest, vfBest, adversarialConfigs, historyRows] = constraintGenerationRobust(
%       modelRef, arm, curveOpts, stressLimit, dispLimit,
%       paramBounds, p0, deltaDeg, cgOpts)
%
%   Minimizes volume fraction subject to stress/displacement constraints holding
%   for all joint-angle configurations in a discrete set.  The adversarial set
%   is discovered incrementally (constraint generation / cutting-plane method):
%
%     outer: solve min VF s.t. constraints for current active config set
%            (via solveCurveParamMinVolume)
%     inner: search worst-case beta over discrete grid on frame surrogate,
%            verify on 3D solid  (via findAdversarialBeta)
%     if new violated config found → add to active set, repeat
%     else → converged (current solution robust over entire discrete space)
%
%   Inputs
%     modelRef    - ManipulatorModel3D reference model (first nominal config)
%     arm         - armModelDefaults struct
%     curveOpts   - options for buildCurveLinkedDensity / evaluateLinkedDensityMetrics
%     stressLimit - stress upper bound [Pa]
%     dispLimit   - tip displacement upper bound [m] (Inf = inactive)
%     paramBounds - struct from defaultCurveParamBounds
%     p0          - initial curve parameter struct
%     deltaDeg    - angular step for discrete beta-grid [deg], e.g. 90 or 45
%     cgOpts      - optional struct:
%                     maxCGIter     max constraint-generation outer iterations (default 10)
%                     topK          candidates per criterion in findAdversarialBeta (default 3)
%                     propLevel     section-prop estimation level 0|1|2 (default 2)
%                     penal         SIMP penalty exponent (default 3)
%                     innerSolverOpts  struct passed to solveCurveParamMinVolume
%                     verbose       print CG progress (default true)
%
%   Outputs
%     pBest            - optimal curve parameter struct
%     vfBest           - optimal volume fraction
%     adversarialConfigs - struct array of all adversarial betas found (with
%                          stressRatio, dispRatio, violated fields)
%     historyRows      - struct array from solveCurveParamMinVolume (last solve)

    if nargin < 9 || isempty(cgOpts)
        cgOpts = struct();
    end
    maxCGIter  = localOpt(cgOpts, 'maxCGIter',  10);
    topK       = localOpt(cgOpts, 'topK',        3);
    propLevel  = localOpt(cgOpts, 'propLevel',   2);
    penal      = localOpt(cgOpts, 'penal',       3);
    verbose    = localOpt(cgOpts, 'verbose',     true);
    innerOpts  = localOpt(cgOpts, 'innerSolverOpts', struct());

    % Default inner solver options if not set.
    innerOpts = applyDefaultOpt(innerOpts, 'maxIter',     80);
    innerOpts = applyDefaultOpt(innerOpts, 'maxFunEvals', 400);
    innerOpts = applyDefaultOpt(innerOpts, 'vfInitial',   1.0);
    innerOpts = applyDefaultOpt(innerOpts, 'vfMin',       0.05);

    % Adversarial-search options forwarded to findAdversarialBeta.
    advOpts.topK       = topK;
    advOpts.penal      = penal;
    advOpts.propLevel  = propLevel;
    advOpts.verbose    = verbose;

    % ------------------------------------------------------------------
    % Seed active set with the nominal straight configuration.
    % ------------------------------------------------------------------
    nominalBeta    = zeros(1, 7);   % first joint fixed at 0, all joints straight
    nominalAnalysis = buildSingleBetaAnalysis(nominalBeta, arm);
    activeAnalyses = {nominalAnalysis};
    activeBetas    = nominalBeta;

    adversarialConfigs = struct([]);
    pCurr  = p0;
    vfCurr = localOpt(innerOpts, 'vfInitial', 1.0);
    historyRows = struct([]);

    if verbose
        fprintf('\n=== Constraint Generation Robust Optimization ===\n');
        fprintf('  stressLimit=%.4e Pa  dispLimit=%.4e m  deltaDeg=%g\n', ...
            stressLimit, dispLimit, deltaDeg);
        fprintf('  maxCGIter=%d  topK=%d  propLevel=%d\n\n', maxCGIter, topK, propLevel);
    end

    % ------------------------------------------------------------------
    % Constraint generation outer loop.
    % ------------------------------------------------------------------
    for cgIter = 1:maxCGIter
        if verbose
            fprintf('--- CG iteration %d/%d  (active configs: %d) ---\n', ...
                cgIter, maxCGIter, numel(activeAnalyses));
        end

        % ---- Inner solve: minimize VF over active config set. --------
        innerOpts.vfInitial = vfCurr;   % warm-start from previous VF
        [pCurr, vfCurr, historyRows] = solveCurveParamMinVolume( ...
            activeAnalyses, modelRef, curveOpts, stressLimit, dispLimit, ...
            paramBounds, pCurr, innerOpts);

        if verbose
            fprintf('  -> VF = %.4f\n', vfCurr);
        end

        % ---- Adversarial search: find worst beta for current design. --
        candidates = findAdversarialBeta(pCurr, vfCurr, arm, modelRef, curveOpts, ...
            stressLimit, dispLimit, deltaDeg, advOpts);

        % Accumulate all candidates in history.
        adversarialConfigs = [adversarialConfigs; candidates(:)]; %#ok<AGROW>

        % Find violated candidates not yet in active set.
        newBetas = selectNewViolatedBetas(candidates, activeBetas);

        if isempty(newBetas)
            if verbose
                fprintf('  No new violated configs found — robust optimum reached.\n\n');
            end
            break;
        end

        % Add new violated configs to active set.
        for ni = 1:size(newBetas, 1)
            betaVec = newBetas(ni, :);
            newAna  = buildSingleBetaAnalysis(betaVec, arm);
            activeAnalyses{end+1} = newAna; %#ok<AGROW>
            activeBetas = [activeBetas; betaVec]; %#ok<AGROW>
            if verbose
                fprintf('  + added config: beta=[%s]\n', num2str(betaVec, '%6.1f '));
            end
        end
    end

    pBest  = pCurr;
    vfBest = vfCurr;

    if verbose
        fprintf('=== Done. Final VF = %.4f, active configs = %d ===\n\n', ...
            vfBest, numel(activeAnalyses));
    end
end

% -------------------------------------------------------------------------
function newBetas = selectNewViolatedBetas(candidates, activeBetas)
% Return rows of violated candidates not already in activeBetas.
% Two beta vectors are considered equal if they match to within 0.1 deg.

    newBetas = zeros(0, size(activeBetas, 2));
    for ci = 1:numel(candidates)
        if ~candidates(ci).violated
            continue;
        end
        betaVec = candidates(ci).beta;
        if ~isAlreadyActive(betaVec, activeBetas)
            newBetas = [newBetas; betaVec]; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
function tf = isAlreadyActive(betaVec, activeBetas)
    if isempty(activeBetas)
        tf = false;
        return;
    end
    diffs = max(abs(activeBetas - betaVec), [], 2);
    tf = any(diffs < 0.1);
end

% -------------------------------------------------------------------------
function analysis = buildSingleBetaAnalysis(betaVec, arm)
    model = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
        arm.res, arm.res_th, arm.alpha, betaVec, arm.ShapeFn, ...
        false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
    analysis = model.analysis;
end

% -------------------------------------------------------------------------
function v = localOpt(s, name, default)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

% -------------------------------------------------------------------------
function s = applyDefaultOpt(s, name, default)
    if ~isfield(s, name)
        s.(name) = default;
    end
end

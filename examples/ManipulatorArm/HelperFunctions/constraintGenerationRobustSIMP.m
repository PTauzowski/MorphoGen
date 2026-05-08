function [optResult, adversarialConfigs, cgHistory] = constraintGenerationRobustSIMP( ...
    baseAnalyses, modelRef, arm, simpOpts, stressLimit, dispLimit, ...
    x0, xmin, xmax, deltaDeg, cgOpts)
% constraintGenerationRobustSIMP  Constraint-generation robust SIMP optimization.
%
%   Wraps solveSIMPVolumeStressMMA in a cutting-plane outer loop that adaptively
%   discovers and adds adversarial beta-configurations as additional load cases.
%
%   Inputs
%     baseAnalyses - {nBase x 1} cell of analyses from armLoadConfigs (always active)
%     modelRef     - ManipulatorModel3D reference for section property estimation
%     arm          - armModelDefaults struct
%     simpOpts     - options struct for solveSIMPVolumeStressMMA (must contain penal,
%                    expand, pullback, configNames, configLabels, etc.)
%     stressLimit  - absolute stress upper bound [Pa] for adversarial violation check
%     dispLimit    - absolute tip-displacement upper bound [m] (Inf = inactive)
%     x0           - initial design vector (optimisation variable space)
%     xmin, xmax   - design variable bounds
%     deltaDeg     - beta-grid angular step [deg] for adversarial search
%     cgOpts       - struct with fields:
%                      maxCGIter  outer CG iterations (default 5)
%                      topK       candidates per criterion in adversarial search (default 3)
%                      propLevel  frame section-prop level (default 1)
%                      penal      SIMP penalty forwarded to findAdversarialBetaSIMP (default simpOpts.penal)
%                      verbose    print progress (default true)
%
%   Outputs
%     optResult          - final solveSIMPVolumeStressMMA result struct (same fields as direct call)
%     adversarialConfigs - struct array of all 3D-FEM-verified adversarial candidates (all CG iters)
%     cgHistory          - struct array, one row per CG iteration

    if nargin < 11 || isempty(cgOpts)
        cgOpts = struct();
    end
    maxCGIter = localOpt(cgOpts, 'maxCGIter', 5);
    topK      = localOpt(cgOpts, 'topK',      3);
    propLevel = localOpt(cgOpts, 'propLevel', 1);
    penal     = localOpt(cgOpts, 'penal',     simpOpts.penal);
    verbose   = localOpt(cgOpts, 'verbose',   true);

    expandFn = localOpt(simpOpts, 'expand', @(z) z);

    % Active set: starts with the caller-supplied base configs.
    % Adversarial betas added by CG are tracked separately to avoid re-adding them.
    activeAnalyses = baseAnalyses(:);
    activeBetas    = zeros(0, 7);   % [nAdv x 7], only CG-discovered betas

    zWarm              = x0(:);
    optResult          = [];
    adversarialConfigs = struct([]);
    cgHistory          = struct([]);

    searchOpts.topK      = topK;
    searchOpts.penal     = penal;
    searchOpts.propLevel = propLevel;
    searchOpts.nJoints   = 7;
    searchOpts.verbose   = verbose;

    for cgIter = 1:maxCGIter
        nActive = numel(activeAnalyses);
        fprintf('\n=== CG iteration %d / %d  (active configs = %d) ===\n', ...
            cgIter, maxCGIter, nActive);

        % Build updated opts with config names matching the current active set.
        currentSimpOpts              = simpOpts;
        currentSimpOpts.configNames  = buildActiveConfigNames(simpOpts, activeBetas);
        currentSimpOpts.configLabels = currentSimpOpts.configNames;

        % Inner SIMP solve (warm-started from previous CG iteration).
        optResult = solveSIMPVolumeStressMMA(activeAnalyses, zWarm, xmin, xmax, currentSimpOpts);
        zWarm = optResult.zFinal;

        vfCurr = optResult.finalVolumeFraction;
        fprintf('  -> VF = %.4f\n', vfCurr);

        % Full-arm physical density for adversarial search.
        rhoFull = expandFn(optResult.zPhysicalFinal);

        % Adversarial beta search.
        candidates = findAdversarialBetaSIMP(rhoFull, arm, modelRef, ...
            stressLimit, dispLimit, deltaDeg, searchOpts);

        % Accumulate all verified candidates across CG iterations.
        for ci = 1:numel(candidates)
            adversarialConfigs = [adversarialConfigs; candidates(ci)]; %#ok<AGROW>
        end

        % Filter to violated betas not already in the active set.
        newBetas = selectNewViolatedBetas(candidates, activeBetas);

        row.cgIter    = cgIter;
        row.nActive   = nActive;
        row.vf        = vfCurr;
        row.nNewBetas = size(newBetas, 1);
        cgHistory = [cgHistory; row]; %#ok<AGROW>

        if isempty(newBetas)
            fprintf('  No new violated configs found — robust optimum reached.\n\n');
            break;
        end

        fprintf('  Adding %d new adversarial config(s) to active set.\n', size(newBetas, 1));
        for bi = 1:size(newBetas, 1)
            betaVec = newBetas(bi, :);
            mdl = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
                arm.res, arm.res_th, arm.alpha, betaVec, arm.ShapeFn, ...
                false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
            activeAnalyses{end+1} = mdl.analysis; %#ok<AGROW>
            activeBetas = [activeBetas; betaVec];  %#ok<AGROW>
        end
    end

    fprintf('\n=== Done. Final VF = %.4f, active configs = %d ===\n', ...
        optResult.finalVolumeFraction, numel(activeAnalyses));
end

% -------------------------------------------------------------------------
function names = buildActiveConfigNames(simpOpts, adversarialBetas)
    if isfield(simpOpts, 'configNames') && ~isempty(simpOpts.configNames)
        baseNames = string(simpOpts.configNames(:));
    else
        baseNames = strings(0, 1);
    end
    nAdv = size(adversarialBetas, 1);
    advNames = strings(nAdv, 1);
    for bi = 1:nAdv
        advNames(bi) = sprintf('adv%02d', bi);
    end
    names = [baseNames; advNames];
end

% -------------------------------------------------------------------------
function newBetas = selectNewViolatedBetas(candidates, activeBetas)
    newBetas = zeros(0, 7);
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
function v = localOpt(s, name, default)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

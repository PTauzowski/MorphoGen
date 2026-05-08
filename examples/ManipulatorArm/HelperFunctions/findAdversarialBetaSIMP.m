function betaCandidates = findAdversarialBetaSIMP(rhoFull, arm, modelRef, ...
    stressLimit, dispLimit, deltaDeg, opts)
% findAdversarialBetaSIMP  Discover adversarial joint-angle configurations for a SIMP density field.
%
%   betaCandidates = findAdversarialBetaSIMP(rhoFull, arm, modelRef,
%       stressLimit, dispLimit, deltaDeg, opts)
%
%   Identical adversarial pipeline to findAdversarialBeta but takes the full-arm
%   density vector directly instead of curve parameters. Suitable for SIMP results.
%
%   Inputs
%     rhoFull     - [nElems x 1] full-arm physical density (SIMP result)
%     arm         - armModelDefaults struct
%     modelRef    - ManipulatorModel3D reference model (used for H and section props)
%     stressLimit - absolute stress upper bound [Pa]
%     dispLimit   - absolute tip-displacement upper bound [m] (Inf = inactive)
%     deltaDeg    - beta-grid angular step [deg]
%     opts        - optional struct:
%                     topK       candidates per criterion (default 3)
%                     penal      SIMP penalty for section-prop estimation (default 3)
%                     propLevel  0|1|2 for estimateFrameSectionPropsFromDensity (default 1)
%                     nJoints    number of beta joints (default 7)
%                     verbose    print progress (default true)
%
%   Output
%     betaCandidates - struct array (sorted by stressRatio descending), fields:
%                       beta, maxHM, tipUz, stressRatio, dispRatio, violated

    if nargin < 7 || isempty(opts)
        opts = struct();
    end
    topK      = localOpt(opts, 'topK',      3);
    penal     = localOpt(opts, 'penal',     3);
    Pz        = localOpt(opts, 'Pz',        arm.Pz);
    propLevel = localOpt(opts, 'propLevel', 1);
    nJoints   = localOpt(opts, 'nJoints',   7);
    verbose   = localOpt(opts, 'verbose',   true);

    rhoFull = rhoFull(:);

    % ------------------------------------------------------------------
    % 1. Extract reference half-segment densities for section prop estimation.
    % ------------------------------------------------------------------
    H      = modelRef.halfSegmentNelems;
    rhoRef = rhoFull(1:H);

    % ------------------------------------------------------------------
    % 2. Density-aware frame section properties.
    % ------------------------------------------------------------------
    if verbose
        fprintf('[findAdversarialBetaSIMP] Estimating frame section props (Level %d)...\n', propLevel);
    end
    props = estimateFrameSectionPropsFromDensity(rhoRef, modelRef, penal, propLevel);

    % ------------------------------------------------------------------
    % 3. Enumerate discrete beta-grid on fast frame model.
    % ------------------------------------------------------------------
    B = buildBetaGrid(nJoints, deltaDeg);
    if verbose
        fprintf('[findAdversarialBetaSIMP] Enumerating %d beta configs on frame...\n', size(B, 1));
    end
    rankSets = enumerateBetaOnFrame(arm, props, B, Pz);

    % ------------------------------------------------------------------
    % 4. Union of top-k indices from all ranking criteria.
    % ------------------------------------------------------------------
    rankFields = {'byAxial','byShearY','byShearZ','byTorsion', ...
                  'byBendingY','byBendingZ','byTipDisp','byComposite'};
    candidateIdx = [];
    for fi = 1:numel(rankFields)
        ranking = rankSets.(rankFields{fi});
        take    = min(topK, numel(ranking));
        candidateIdx = union(candidateIdx, ranking(1:take));
    end
    candidateIdx = candidateIdx(:);
    if verbose
        fprintf('[findAdversarialBetaSIMP] %d unique frame-level candidates selected.\n', ...
            numel(candidateIdx));
    end

    % ------------------------------------------------------------------
    % 5. Verify candidates on full 3D solid FEM.
    % ------------------------------------------------------------------
    nCand = numel(candidateIdx);
    betaCandidates = struct( ...
        'beta',        cell(nCand, 1), ...
        'maxHM',       cell(nCand, 1), ...
        'tipUz',       cell(nCand, 1), ...
        'stressRatio', cell(nCand, 1), ...
        'dispRatio',   cell(nCand, 1), ...
        'violated',    cell(nCand, 1));

    dispRef        = max(abs(dispLimit), eps);
    candidateBetas = B(candidateIdx, :);

    sRatioVec   = zeros(nCand, 1);
    dRatioVec   = zeros(nCand, 1);
    maxHMVec    = zeros(nCand, 1);
    tipUzVec    = zeros(nCand, 1);
    violatedVec = false(nCand, 1);

    evalOpts       = struct();
    evalOpts.penal = penal;

    parfor ci = 1:nCand
        betaVec = candidateBetas(ci, :);
        mdl = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
            arm.res, arm.res_th, arm.alpha, betaVec, arm.ShapeFn, ...
            false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
        [metrics, ~] = evaluateLinkedDensityMetrics({mdl.analysis}, rhoFull, evalOpts);

        sR = metrics.maxHM(1) / stressLimit;
        dR = abs(metrics.tipUz(1)) / dispRef;
        sRatioVec(ci)   = sR;
        dRatioVec(ci)   = dR;
        maxHMVec(ci)    = metrics.maxHM(1);
        tipUzVec(ci)    = metrics.tipUz(1);
        violatedVec(ci) = (sR > 1) || (isfinite(dispLimit) && dR > 1);
    end

    for ci = 1:nCand
        betaCandidates(ci).beta        = candidateBetas(ci, :);
        betaCandidates(ci).maxHM       = maxHMVec(ci);
        betaCandidates(ci).tipUz       = tipUzVec(ci);
        betaCandidates(ci).stressRatio = sRatioVec(ci);
        betaCandidates(ci).dispRatio   = dRatioVec(ci);
        betaCandidates(ci).violated    = violatedVec(ci);
    end

    if verbose
        for ci = 1:nCand
            if violatedVec(ci), tag = 'VIOLATED'; else, tag = 'ok'; end
            fprintf('  cand %2d/%2d: beta=[%s]  s=%.3f  d=%.3f  [%s]\n', ...
                ci, nCand, num2str(candidateBetas(ci,:), '%6.1f '), ...
                sRatioVec(ci), dRatioVec(ci), tag);
        end
    end

    % Sort by stress ratio descending (worst first).
    sRatios = [betaCandidates.stressRatio];
    [~, ord] = sort(sRatios, 'descend');
    betaCandidates = betaCandidates(ord);
end

% -------------------------------------------------------------------------
function v = localOpt(s, name, default)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

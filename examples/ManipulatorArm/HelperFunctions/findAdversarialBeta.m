function betaCandidates = findAdversarialBeta(p, vf, arm, modelRef, curveOpts, ...
    stressLimit, dispLimit, deltaDeg, opts)
% findAdversarialBeta  Discover adversarial joint-angle configurations.
%
%   betaCandidates = findAdversarialBeta(p, vf, arm, modelRef, curveOpts,
%       stressLimit, dispLimit, deltaDeg, opts)
%
%   Given a fixed density field (curve params p at volume fraction vf), finds
%   the worst-case joint-angle vectors beta from a discrete grid by:
%     1. Estimating density-aware frame section properties (Level-2 surrogate).
%     2. Enumerating the full discrete beta-grid on the fast frame model.
%     3. Taking the union of top-k candidates from 8 ranking criteria.
%     4. Verifying every candidate on the full 3D solid FEM.
%   Returns candidates sorted by stress ratio descending (worst first).
%
%   Inputs
%     p           - curve parameter struct (passed to buildCurveLinkedDensity)
%     vf          - volume fraction scalar (overrides curveOpts.VolFrac)
%     arm         - armModelDefaults struct (E, nu, R, r, h_seg, res, res_th,
%                   alpha, Pz, ShapeFn, constEndRing, constMiddleRing, nCircDiv)
%     modelRef    - ManipulatorModel3D reference model (first config)
%     curveOpts   - options struct for buildCurveLinkedDensity /
%                   evaluateLinkedDensityMetrics
%     stressLimit - stress upper bound [Pa]
%     dispLimit   - tip displacement upper bound [m] (Inf = inactive)
%     deltaDeg    - angular step for beta-grid [deg], e.g. 90 or 45
%     opts        - optional struct:
%                     topK        number of candidates per criterion (default 3)
%                     penal       SIMP penalty for section-prop estimation (default 3)
%                     Pz          frame tip load for enumeration [N] (default arm.Pz)
%                     propLevel   0|1|2 for estimateFrameSectionPropsFromDensity (default 2)
%                     nJoints     number of beta joints (default 7)
%                     verbose     print progress (default true)
%
%   Output
%     betaCandidates - struct array, one entry per verified candidate, fields:
%                       beta         [1 x nJoints] joint angles [deg]
%                       maxHM        max von-Mises stress [Pa]
%                       tipUz        tip displacement [m]
%                       stressRatio  maxHM / stressLimit
%                       dispRatio    |tipUz| / |dispLimit|
%                       violated     true if stressRatio>1 or (finite dispLimit && dispRatio>1)

    if nargin < 9 || isempty(opts)
        opts = struct();
    end
    topK       = localOpt(opts, 'topK',      3);
    penal      = localOpt(opts, 'penal',     3);
    Pz         = localOpt(opts, 'Pz',        arm.Pz);
    propLevel  = localOpt(opts, 'propLevel', 1);
    nJoints    = localOpt(opts, 'nJoints',   7);
    verbose    = localOpt(opts, 'verbose',   true);

    % ------------------------------------------------------------------
    % 1. Build density field for the current (p, vf) point.
    % ------------------------------------------------------------------
    localCurveOpts         = curveOpts;
    localCurveOpts.VolFrac = vf;
    [rhoRef, rhoFull, ~] = buildCurveLinkedDensity(p, modelRef, localCurveOpts);

    % ------------------------------------------------------------------
    % 2. Density-aware frame section properties (reference half-segment only).
    % ------------------------------------------------------------------
    if verbose
        fprintf('[findAdversarialBeta] Estimating frame section props (Level %d)...\n', propLevel);
    end
    props = estimateFrameSectionPropsFromDensity(rhoRef, modelRef, penal, propLevel);

    % ------------------------------------------------------------------
    % 3. Enumerate discrete beta-grid on fast frame model.
    % ------------------------------------------------------------------
    B = buildBetaGrid(nJoints, deltaDeg);   % [N x nJoints]
    if verbose
        fprintf('[findAdversarialBeta] Enumerating %d beta configs on frame...\n', size(B, 1));
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
        fprintf('[findAdversarialBeta] %d unique frame-level candidates selected.\n', ...
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

    dispRef = max(abs(dispLimit), eps);

    % Pre-extract candidate beta rows so parfor has no shared indexing.
    candidateBetas = B(candidateIdx, :);

    sRatioVec   = zeros(nCand, 1);
    dRatioVec   = zeros(nCand, 1);
    maxHMVec    = zeros(nCand, 1);
    tipUzVec    = zeros(nCand, 1);
    violatedVec = false(nCand, 1);

    parfor ci = 1:nCand
        betaVec  = candidateBetas(ci, :);
        analysis = buildSingleBetaAnalysis(betaVec, arm);
        [metrics, ~] = evaluateLinkedDensityMetrics({analysis}, rhoFull, localCurveOpts);

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
function analysis = buildSingleBetaAnalysis(betaVec, arm)
% Build a LinearElasticityWeighted for a single beta configuration.
% Mirrors the model construction in testCurveParamSegmentOptimization.

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

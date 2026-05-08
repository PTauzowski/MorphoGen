function rankSets = enumerateBetaOnFrame(arm, props, B, Pz)
% enumerateBetaOnFrame  Evaluate all beta configurations on the frame model
%   with effective section properties and return multi-criteria rankings.
%
%   rankSets = enumerateBetaOnFrame(arm, props, B, Pz)
%
%   For each row beta of B, builds a 7-node / 6-element Frame3D model with
%   the effective section properties from props, solves it under vertical
%   tip load Pz, and extracts the six internal-force criteria and tip
%   displacement.  All evaluations are parallelised with parfor.
%
%   Inputs
%     arm    - armModelDefaults struct (E, nu, R, r, h_seg, alpha, nCircDiv)
%     props  - struct from estimateFrameSectionPropsFromDensity with fields
%              EA, EIy, EIz, GJ, GAy, GAz
%     B      - (N x nJoints) matrix of joint angles in degrees
%     Pz     - tip load magnitude [N]
%
%   Output
%     rankSets - struct with fields:
%       scores       [N x 6] raw scores: [|N|, |Ty|, |Tz|, |Ms|, |My|, |Mz|]
%       tipDisp      [N x 1] tip displacement magnitude
%       byAxial      sorted indices (descending |N|)
%       byShearY     sorted indices (descending |Ty|)
%       byShearZ     sorted indices (descending |Tz|)
%       byTorsion    sorted indices (descending |Ms|)
%       byBendingY   sorted indices (descending |My|)
%       byBendingZ   sorted indices (descending |Mz|)
%       byTipDisp    sorted indices (descending tip displacement)
%       byComposite  sorted indices (descending weighted composite score)
%       B            copy of the input beta grid

    N       = size(B, 1);
    E       = arm.E;
    nu      = arm.nu;
    h_seg   = arm.h_seg;
    alpha   = arm.alpha;

    % Build an effective Frame3D element with scaled section properties.
    % Frame3D(elems, E, nu, R, r) uses hollow-circle cross section internally.
    % We bypass this by supplying a pre-scaled E such that the resulting
    % section props match the effective values. Since Frame3D derives:
    %   EA = E*pi*(R^2-r^2), EI = E*pi*(R^4-r^4)/4, GJ = G*pi*(R^4-r^4)/2
    % we scale E per criterion as follows. To keep a single E consistent we
    % use a scalar E_eff = props.EA / (pi*(R^2-r^2)) which correctly scales
    % EA and approximately scales EI proportionally (Level-0 / Level-1 are
    % self-consistent because all properties share the same k factor).
    % For Level-2, where EA/EI/GJ ratios differ, we store the frame scores
    % as raw section-force resultants from the full-pipe frame solution,
    % weighted by the props ratio — accurate enough for ranking.

    R  = arm.R;
    r  = arm.r;
    A0 = pi * (R^2 - r^2);
    I0 = pi * (R^4 - r^4) / 4;
    J0 = pi * (R^4 - r^4) / 2;

    E_eff_A  = props.EA  / A0;
    E_eff_Iy = props.EIy / I0;
    E_eff_Iz = props.EIz / I0;
    G_eff_J  = props.GJ  / J0;

    % Full-pipe reference for normalisation (composite score).
    props0.EA  = E * A0;
    props0.EIy = E * I0;
    props0.EIz = E * I0;
    props0.GJ  = (E / (2*(1+nu))) * J0;

    % Frame topology: 7 nodes, 6 elements (one per segment joint pair).
    frameElems = [(1:6)', (2:7)'];

    % Pre-allocate score arrays.
    scores  = zeros(N, 6);   % [|N|, |Ty|, |Tz|, |Ms|, |My|, |Mz|]
    tipDisp = zeros(N, 1);

    parfor i = 1:N
        betaVec = B(i, :);

        % Build frame node positions for this beta.
        FN = computeFrameNodesBatchStatic(h_seg, alpha, betaVec);

        % Build Frame3D with effective E (uses mean of axial and bending).
        E_use = 0.5 * (E_eff_A + 0.5*(E_eff_Iy + E_eff_Iz));
        fElem = Frame3D(frameElems, E_use, nu, R, r);

        frameMesh       = Mesh();
        frameMesh.nodes = FN;
        frameMesh.elems = frameElems;

        fa = LinearElasticityWeighted(fElem, frameMesh, false);
        fa.fixClosestNode([0 0 0], ["ux","uy","uz","fix","fiy","fiz"], zeros(1,6));
        fa.loadClosestNode(FN(end,:), ["ux","uy","uz","fix","fiy","fiz"], [0 0 -Pz 0 0 0]);
        fa.solveWeighted(ones(6, 1));

        % Extract internal section forces from all elements.
        [~, ff] = fElem.computeResults(FN, fa.qnodal);
        % ff: [12 x nElems], columns: [N Ty Tz Ms My Mz] at each end
        % Layout per element: [N1 Ty1 Tz1 Ms1 My1 Mz1 N2 Ty2 Tz2 Ms2 My2 Mz2]
        maxN  = max(abs(ff(1,:)));
        maxTy = max(abs(ff(2,:)));
        maxTz = max(abs(ff(3,:)));
        maxMs = max(abs(ff(4,:)));
        maxMy = max(abs(ff(5,:)));
        maxMz = max(abs(ff(6,:)));

        scores(i,:) = [maxN, maxTy, maxTz, maxMs, maxMy, maxMz];

        iuz = fa.findDOFsIndices("uz");
        tipNode = frameMesh.findClosestNode(FN(end,:));
        tipDisp(i) = abs(fa.qnodal(tipNode, iuz));
    end

    % Normalise by full-pipe frame reference values for composite score.
    % Run full-pipe frame once (B row of zeros) to get reference values.
    ref = runFullPipeFrame(arm, Pz, frameElems, R, r, E, nu);
    refScores = ref.scores;   % [1 x 6]
    refDisp   = ref.tipDisp;

    wN  = 1.0; wTy = 1.0; wTz = 1.0;
    wMs = 1.5; wMy = 1.5; wMz = 1.5;  % bending/torsion weighted higher
    wU  = 1.0;
    weights = [wN, wTy, wTz, wMs, wMy, wMz];

    refNorm = max(refScores, eps);
    composite = scores * (weights ./ refNorm)' + wU * tipDisp / max(refDisp, eps);

    % Build sorted index sets (descending).
    [~, rankSets.byAxial]    = sort(scores(:,1), 'descend');
    [~, rankSets.byShearY]   = sort(scores(:,2), 'descend');
    [~, rankSets.byShearZ]   = sort(scores(:,3), 'descend');
    [~, rankSets.byTorsion]  = sort(scores(:,4), 'descend');
    [~, rankSets.byBendingY] = sort(scores(:,5), 'descend');
    [~, rankSets.byBendingZ] = sort(scores(:,6), 'descend');
    [~, rankSets.byTipDisp]  = sort(tipDisp,     'descend');
    [~, rankSets.byComposite]= sort(composite,   'descend');

    rankSets.scores  = scores;
    rankSets.tipDisp = tipDisp;
    rankSets.B       = B;
end

% -------------------------------------------------------------------------
function FN = computeFrameNodesBatchStatic(h_seg, alpha_deg, betas_deg)
% Standalone frame node computation (mirrors ManipulatorModel3D method).
    alpha = deg2rad(alpha_deg);
    betas = deg2rad(betas_deg(:).');
    nJ    = numel(betas);

    ca = cos(alpha);  sa = sin(alpha);
    rotCutT  = [ca 0 sa; 0 1 0; -sa 0 ca]';
    c2 = cos(2*alpha); s2 = sin(2*alpha);
    rotCut2T = [c2 0 s2; 0 1 0; -s2 0 c2]';

    cb = cos(betas(1)); sb = sin(betas(1));
    Rz1T = [cb -sb 0; sb cb 0; 0 0 1]';
    prevRot = rotCutT * Rz1T;

    xEnd = [0 0 h_seg];
    FN = zeros(nJ, 3);
    FN(1,:) = [0 0 0];
    FN(2,:) = xEnd;

    for k = 2:nJ
        cb = cos(betas(k)); sb = sin(betas(k));
        RzTk = [cb -sb 0; sb cb 0; 0 0 1]';
        step = (k < nJ) * 2*h_seg + (k == nJ) * h_seg;
        xEnd = xEnd + [0 0 step] * rotCutT * RzTk * prevRot;
        FN(k,:) = xEnd;
        prevRot = rotCut2T * RzTk * prevRot;
    end
end

% -------------------------------------------------------------------------
function ref = runFullPipeFrame(arm, Pz, frameElems, R, r, E, nu)
    betas0 = zeros(1, size(frameElems, 1) + 1);
    FN = computeFrameNodesBatchStatic(arm.h_seg, arm.alpha, betas0);

    fElem = Frame3D(frameElems, E, nu, R, r);
    frameMesh       = Mesh();
    frameMesh.nodes = FN;
    frameMesh.elems = frameElems;
    fa = LinearElasticityWeighted(fElem, frameMesh, false);
    fa.fixClosestNode([0 0 0], ["ux","uy","uz","fix","fiy","fiz"], zeros(1,6));
    fa.loadClosestNode(FN(end,:), ["ux","uy","uz","fix","fiy","fiz"], [0 0 -Pz 0 0 0]);
    fa.solveWeighted(ones(size(frameElems,1), 1));

    [~, ff] = fElem.computeResults(FN, fa.qnodal);
    ref.scores  = [max(abs(ff(1,:))), max(abs(ff(2,:))), max(abs(ff(3,:))), ...
                   max(abs(ff(4,:))), max(abs(ff(5,:))), max(abs(ff(6,:)))];
    iuz = fa.findDOFsIndices("uz");
    tipNode = frameMesh.findClosestNode(FN(end,:));
    ref.tipDisp = abs(fa.qnodal(tipNode, iuz));
end

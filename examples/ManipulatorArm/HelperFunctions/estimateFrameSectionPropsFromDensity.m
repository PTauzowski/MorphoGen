function props = estimateFrameSectionPropsFromDensity(rhoRef, model, penal, level)
% estimateFrameSectionPropsFromDensity  Surrogate effective beam section
%   properties for a density-penalized Arm-Z reference module.
%
%   props = estimateFrameSectionPropsFromDensity(rhoRef, model, penal, level)
%
%   This function is a SURROGATE, not a true homogenization. The three
%   levels trade accuracy for computational cost:
%
%   Level 0 — scalar scaling
%     All section properties scaled uniformly by mean(rho^p).
%     Use for coarse screening or early CG iterations.
%
%   Level 1 — section-moment projection  (engineering surrogate)
%     EA, EIy, EIz, GJ estimated from cross-section integrals of rho^p.
%     Captures relative variation between bending, torsion, and axial
%     stiffness, but assumes a straight prismatic cross-section perpendicular
%     to the local beam axis. Valid as a surrogate; do not call it
%     homogenization in a paper.
%
%   Level 2 — numerical beam-equivalent stiffness  (default for paper runs)
%     Six unit boundary-value problems on the penalized 3D reference module
%     with Saint-Venant style BCs (resultant coupling, free warping).
%     Loads and reactions expressed in local beam axis frame.
%     This is the defensible level for publication.
%
%   Inputs
%     rhoRef  - [H x 1] reference half-segment element densities
%     model   - ManipulatorModel3D reference model
%     penal   - SIMP penalty exponent
%     level   - 0, 1, or 2  (default: 2)
%
%   Output
%     props   - struct with fields:
%                 EA, EIy, EIz, GJ, GAy, GAz  (effective, local frame)
%                 EA0, EIy0, EIz0, GJ0         (full-pipe reference values)
%                 level                         (which level was used)

    if nargin < 4 || isempty(level)
        level = 2;
    end

    E  = model.fe.mat.E;
    nu = model.fe.mat.nu;
    G  = E / (2 * (1 + nu));
    R  = model.R;
    r  = model.r;

    % Cowper (1966) shear correction factor for hollow circular cross-section.
    m_r = r / R;
    kappa = 6*(1+nu)*(1+m_r^2)^2 / ...
            ((7+6*nu)*(1+m_r^2)^2 + (20+12*nu)*m_r^2);

    % Full-pipe (solid annular) reference section properties.
    % GAy0 / GAz0 include kappa so they represent the Timoshenko effective
    % shear stiffness used directly in Frame3DSectionProps.
    props.EA0  = E * pi * (R^2 - r^2);
    props.EIy0 = E * pi * (R^4 - r^4) / 4;
    props.EIz0 = props.EIy0;
    props.GJ0  = G * pi * (R^4 - r^4) / 2;
    props.GAy0 = kappa * G * pi * (R^2 - r^2);
    props.GAz0 = props.GAy0;
    props.level = level;

    switch level
        case 0
            props = scalarScaling(props, rhoRef, penal);
        case 1
            props = sectionMomentProjection(props, rhoRef, model, penal, E, G, kappa);
        case 2
            props = numericalBeamEquivalent(props, rhoRef, model, penal);
        otherwise
            error('estimateFrameSectionPropsFromDensity: level must be 0, 1, or 2.');
    end
end

% =========================================================================
function props = scalarScaling(props, rhoRef, penal)
    k = mean(rhoRef(:) .^ penal);
    props.EA  = k * props.EA0;
    props.EIy = k * props.EIy0;
    props.EIz = k * props.EIz0;
    props.GJ  = k * props.GJ0;
    props.GAy = k * props.GAy0;
    props.GAz = props.GAy;
end

% =========================================================================
function props = sectionMomentProjection(props, rhoRef, model, penal, E, G, kappa)
% Section-moment projection in the local beam frame of the reference segment.
% The local beam axis x is along the first segment (from node 1 to node 2
% of the frame model). y and z are the cross-section axes.

    H     = model.halfSegmentNelems;
    elems = model.mesh.elems(1:H, :);
    nodes = model.mesh.nodes;

    % Element centroids.
    centroids = squeeze(mean(reshape(nodes(elems', :), ...
        size(elems, 2), H, 3), 1));   % [H x 3]

    % Local beam axis x: direction of first frame segment.
    xA = model.frameNodes(1, :);
    xB = model.frameNodes(2, :);
    ex = (xB - xA) / norm(xB - xA);

    % Build orthonormal local frame (ex, ey, ez).
    % Choose ey perpendicular to ex in the plane containing global z if possible.
    globalZ = [0 0 1];
    if abs(dot(ex, globalZ)) > 0.99
        globalZ = [0 1 0];
    end
    ez = cross(ex, globalZ);
    ez = ez / norm(ez);
    ey = cross(ez, ex);
    ey = ey / norm(ey);

    % Project centroids onto local (ey, ez) cross-section axes.
    c0 = mean(centroids, 1);  % approximate centroid of cross-section
    dc = centroids - c0;
    y_local = dc * ey';   % [H x 1]
    z_local = dc * ez';   % [H x 1]

    % Approximate element cross-sectional area from mesh volume / segment length.
    L = norm(xB - xA);
    totalArea = pi * (model.R^2 - model.r^2);
    Ai = (totalArea / H) * ones(H, 1);   % uniform element area share

    w = rhoRef(:) .^ penal;  % penalized density weights [H x 1]

    props.EA  = E     * sum(w .* Ai);
    props.EIy = E     * sum(w .* z_local.^2 .* Ai);
    props.EIz = E     * sum(w .* y_local.^2 .* Ai);
    props.GJ  = G     * sum(w .* (y_local.^2 + z_local.^2) .* Ai);
    props.GAy = kappa * G * sum(w .* Ai);
    props.GAz = props.GAy;
    props.L   = L;
end

% =========================================================================
function props = numericalBeamEquivalent(props, rhoRef, model, penal)
% Numerical beam-equivalent stiffness via 6 unit BVPs on the penalized 3D
% reference half-segment with Saint-Venant boundary conditions.
%
% Only the reference half-segment (elements 1:H) is penalized.  All other
% arm elements are assigned zero density so they contribute no stiffness
% and do not contaminate the section property estimate through shared nodes.
%
% BVP setup per load case:
%   Root face: all DOFs fixed (identified from model.analysis.supports).
%   Tip face:  distributed unit resultant; nodes free to warp.
%   Non-half-segment DOFs excluded implicitly (freeDofs ⊂ halfSegDofs).
% Six BVPs solved simultaneously (one factorisation, 6 RHS).
%
% Moment loads (BVPs 4-6) are normalised after assembly so the resultant
% moment about the respective local axis equals exactly 1 N·m.
%
% Stiffness inversion:
%   raw EA, GJ, EIy, EIz are recovered from generalized tip translation or
%   least-squares fitted tip-section rotation.  raw GAy and GAz subtract the
%   Euler-Bernoulli bending contribution from transverse-force BVPs.
%   The raw values are calibrated by the same BVP applied to rho=1, then
%   scaled by the analytical annular full-pipe properties used by Frame3D.

    H     = model.halfSegmentNelems;
    elems = model.mesh.elems(1:H, :);

    % Local beam axis direction from first frame segment.
    xA = model.frameNodes(1, :);
    xB = model.frameNodes(2, :);
    ex = (xB - xA) / norm(xB - xA);

    globalZ = [0 0 1];
    if abs(dot(ex, globalZ)) > 0.99, globalZ = [0 1 0]; end
    ez_loc = cross(ex, globalZ); ez_loc = ez_loc / norm(ez_loc);
    ey_loc = cross(ez_loc, ex);  ey_loc = ey_loc / norm(ey_loc);
    R_loc  = [ex; ey_loc; ez_loc];   % [3×3] global→local

    % ------------------------------------------------------------------
    % Half-segment nodes and exact end cross-sections.
    % ------------------------------------------------------------------
    halfSegNodeIds = unique(elems(:));
    L              = norm(xB - xA);
    [rootNodeIds, tipNodeIds] = referenceHalfSegmentSections(model);
    rootNodeIds = intersect(int32(rootNodeIds(:)), int32(halfSegNodeIds(:)));
    tipNodeIds  = intersect(int32(tipNodeIds(:)),  int32(halfSegNodeIds(:)));
    assert(~isempty(rootNodeIds) && ~isempty(tipNodeIds), ...
        'estimateFrameSectionPropsFromDensity: could not identify reference half-segment end sections.');

    raw = solveBeamEquivalentRaw(model, rhoRef, penal, halfSegNodeIds, ...
        rootNodeIds, tipNodeIds, ex, ey_loc, ez_loc, R_loc, L);
    fullRaw = solveBeamEquivalentRaw(model, ones(H, 1), 1.0, halfSegNodeIds, ...
        rootNodeIds, tipNodeIds, ex, ey_loc, ez_loc, R_loc, L);

    % Calibrate the BVP against its own full-pipe response.  The first
    % half-segment has angled Arm-Z end sections, so the raw BVP is used for
    % relative density sensitivity while the full-pipe scale is anchored to
    % the analytical annular properties consumed by the frame element.
    props.EA  = props.EA0  * raw.EA  / max(fullRaw.EA,  eps);
    props.EIy = props.EIy0 * raw.EIy / max(fullRaw.EIy, eps);
    props.EIz = props.EIz0 * raw.EIz / max(fullRaw.EIz, eps);
    props.GJ  = props.GJ0  * raw.GJ  / max(fullRaw.GJ,  eps);
    props.GAy = props.GAy0 * raw.GAy / max(fullRaw.GAy, eps);
    props.GAz = props.GAz0 * raw.GAz / max(fullRaw.GAz, eps);
    props.L   = L;
end

% =========================================================================
function raw = solveBeamEquivalentRaw(model, rhoRef, penal, halfSegNodeIds, ...
    rootNodeIds, tipNodeIds, ex, ey_loc, ez_loc, R_loc, L)

    nodes = model.mesh.nodes;

    % DOF layout: 3 DOFs per node, interleaved [ux_1,uy_1,uz_1, ux_2,...].
    nDPNode = 3;
    nDof    = size(nodes, 1) * nDPNode;
    dofOf   = @(ids, d) (ids(:) - 1) * nDPNode + d;   % d=1:ux, 2:uy, 3:uz

    halfSegDofs = unique([dofOf(halfSegNodeIds,1); dofOf(halfSegNodeIds,2); dofOf(halfSegNodeIds,3)]);
    rootDofs    = unique([dofOf(rootNodeIds,1);    dofOf(rootNodeIds,2);    dofOf(rootNodeIds,3)]);
    freeDofs    = setdiff(halfSegDofs, rootDofs);

    xPenal     = rhoRef(:) .^ penal;
    nElemsFull = model.analysis.getTotalElemsNumber();
    xPenalFull = zeros(nElemsFull, 1);
    xPenalFull(1:numel(rhoRef)) = xPenal;

    [I, J, ~] = model.analysis.globalMatrixIndices();
    stiffFn   = model.analysis.weightedStiffnessFunction();
    Kvals     = model.analysis.globalMatrixAggregationWeighted(stiffFn, xPenalFull);
    K         = sparse(I, J, Kvals, nDof, nDof);
    Kff       = K(freeDofs, freeDofs);

    tipCentroid = mean(nodes(tipNodeIds, :), 1);
    nTipNodes   = numel(tipNodeIds);
    F_full      = zeros(nDof, 6);

    for iLoad = 1:6
        for iN = 1:nTipNodes
            nid   = tipNodeIds(iN);
            r_vec = nodes(nid, :) - tipCentroid;
            switch iLoad
                case 1,  fGlob = ex      / nTipNodes;
                case 2,  fGlob = ey_loc  / nTipNodes;
                case 3,  fGlob = ez_loc  / nTipNodes;
                case 4,  fGlob = cross(ex,     r_vec);
                case 5,  fGlob = cross(ey_loc, r_vec);
                case 6,  fGlob = cross(ez_loc, r_vec);
            end
            F_full(dofOf(nid,1), iLoad) = F_full(dofOf(nid,1), iLoad) + fGlob(1);
            F_full(dofOf(nid,2), iLoad) = F_full(dofOf(nid,2), iLoad) + fGlob(2);
            F_full(dofOf(nid,3), iLoad) = F_full(dofOf(nid,3), iLoad) + fGlob(3);
        end
        switch iLoad
            case 4
                F_full(:,4) = normalizeMomentLoad(F_full(:,4), nodes, tipNodeIds, tipCentroid, ex);
            case 5
                F_full(:,5) = normalizeMomentLoad(F_full(:,5), nodes, tipNodeIds, tipCentroid, ey_loc);
            case 6
                F_full(:,6) = normalizeMomentLoad(F_full(:,6), nodes, tipNodeIds, tipCentroid, ez_loc);
        end
    end

    U_free = Kff \ F_full(freeDofs, :);
    U = zeros(nDof, 6);
    U(freeDofs, :) = U_free;

    [tipU_loc, tipTheta_loc] = fitTipSectionRigidMotion( ...
        U, nodes, tipNodeIds, tipCentroid, R_loc);

    raw.EA = L / max(abs(tipU_loc(1,1)), eps);
    raw.GJ = L / max(abs(tipTheta_loc(1,4)), eps);
    raw.EIy = L / max(abs(tipTheta_loc(2,5)), eps);
    raw.EIz = L / max(abs(tipTheta_loc(3,6)), eps);

    u_y_bend = L^3 / (3 * max(raw.EIz, eps));
    raw.GAy  = L / max(abs(tipU_loc(2,2)) - u_y_bend, eps);

    u_z_bend = L^3 / (3 * max(raw.EIy, eps));
    raw.GAz  = L / max(abs(tipU_loc(3,3)) - u_z_bend, eps);
end

% =========================================================================
function [uLoc, thetaLoc] = fitTipSectionRigidMotion(U, nodes, tipNodeIds, tipCentroid, R_loc)
% Least-squares section motion u_i = u0 + theta × (x_i - c).
    tipNodeIds = tipNodeIds(:);
    nTip = numel(tipNodeIds);
    A = zeros(3*nTip, 6);
    for i = 1:nTip
        nid = tipNodeIds(i);
        r = nodes(nid, :) - tipCentroid;
        rows = 3*(i-1) + (1:3);
        A(rows, 1:3) = eye(3);
        A(rows, 4:6) = -skewMatrix(r);
    end

    uLoc = zeros(3, size(U, 2));
    thetaLoc = zeros(3, size(U, 2));
    for loadCase = 1:size(U, 2)
        b = zeros(3*nTip, 1);
        for i = 1:nTip
            nid = tipNodeIds(i);
            dofs = [3*(nid-1)+1, 3*(nid-1)+2, 3*(nid-1)+3];
            b(3*(i-1) + (1:3)) = U(dofs, loadCase);
        end
        q = A \ b;
        uLoc(:, loadCase) = R_loc * q(1:3);
        thetaLoc(:, loadCase) = R_loc * q(4:6);
    end
end

% =========================================================================
function S = skewMatrix(r)
    S = [  0,   -r(3),  r(2); ...
          r(3),  0,    -r(1); ...
         -r(2), r(1),   0   ];
end

% =========================================================================
function [rootNodeIds, tipNodeIds] = referenceHalfSegmentSections(model)
% Exact section node sets captured during ManipulatorModel3D construction.
    rootNodeIds = [];
    tipNodeIds = [];
    if isprop(model, 'couplingSections') && ~isempty(model.couplingSections)
        secs = model.couplingSections;
        for si = 1:numel(secs)
            if secs(si).frameNode == 1
                rootNodeIds = secs(si).nodes;
            elseif secs(si).frameNode == 2 && secs(si).adjacentFrameNode == 1
                tipNodeIds = secs(si).nodes;
            end
        end
    end
    if isempty(rootNodeIds)
        rootNodeIds = find(any(model.analysis.supports, 2));
    end
    if isempty(tipNodeIds)
        H = model.halfSegmentNelems;
        elems = model.mesh.elems(1:H, :);
        nodes = model.mesh.nodes;
        xA = model.frameNodes(1, :);
        xB = model.frameNodes(2, :);
        ex = (xB - xA) / norm(xB - xA);
        halfSegNodeIds = unique(elems(:));
        xProj = nodes(halfSegNodeIds, :) * ex';
        % Tolerance must cover the full angled cut face: 2·R·sin(alpha).
        tol = 2 * model.R * sin(model.alpha) * 1.1 + norm(xB - xA) * 0.01;
        tipNodeIds = halfSegNodeIds(xProj >= max(xProj) - tol);
    end
end

% =========================================================================
function F = normalizeMomentLoad(F, nodes, tipNodeIds, tipCentroid, axis)
% Remove any numerical net force, then scale resultant moment to one.
    tipNodeIds = tipNodeIds(:);
    nTip = numel(tipNodeIds);
    fMean = [0 0 0];
    for i = 1:nTip
        nid = tipNodeIds(i);
        dofs = [3*(nid-1)+1, 3*(nid-1)+2, 3*(nid-1)+3];
        fMean = fMean + F(dofs)';
    end
    fMean = fMean / max(nTip, 1);
    for i = 1:nTip
        nid = tipNodeIds(i);
        dofs = [3*(nid-1)+1, 3*(nid-1)+2, 3*(nid-1)+3];
        F(dofs) = F(dofs) - fMean';
    end
    moment = computeResultantTorque(F, nodes, tipNodeIds, tipCentroid, axis);
    assert(abs(moment) > eps, ...
        'estimateFrameSectionPropsFromDensity: zero resultant moment in section-property BVP.');
    F = F / moment;
end

% =========================================================================
function T = computeResultantTorque(F, nodes, tipNodeIds, tipCentroid, ex)
    T = 0;
    for i = 1:numel(tipNodeIds)
        nid = tipNodeIds(i);
        r_vec = nodes(nid, :) - tipCentroid;
        f_vec = F([3*(nid-1)+1, 3*(nid-1)+2, 3*(nid-1)+3])';
        T = T + dot(cross(r_vec, f_vec), ex);
    end
end

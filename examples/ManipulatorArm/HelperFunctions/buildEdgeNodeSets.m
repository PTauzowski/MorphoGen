function edgeSets = buildEdgeNodeSets(solidNodes, frameNodes, frameElems, R_outer, alpha_rad)
% BUILDEDGENODESETS  Identify solid mesh nodes belonging to each frame-node
% cross-section for kinematic coupling in a multi-segment manipulator arm.
%
% This helper mirrors the class-based selector in ManipulatorModel3D but
% operates on the supplied node cloud only.  To avoid constraining thick
% volumetric bands, the caller should pass boundary nodes here.  The kept
% alpha input is unused and only preserved for API compatibility.
%
% INPUTS
%   solidNodes  - nS x 3  solid boundary-node coordinates  [m]
%   frameNodes  - nF x 3  frame joint positions         [m]
%   frameElems  - nE x 2  frame element connectivity    (1-based node indices)
%   R_outer     - scalar outer tube radius            [m]
%   alpha_rad   - scalar cross-section inclination angle [rad]
%
% OUTPUT
%   edgeSets    - {nF x 1} cell array; edgeSets{k} is a column vector of
%                 solid node indices (1-based) belonging to frame node k.

nFN = size(frameNodes, 1);
boundaryNodeIds = int32((1:size(solidNodes, 1)).');
boundaryCoords = solidNodes;
snapTol = 1.0e-6;
radialTol = R_outer * 1.05 + snapTol;

% Kept for compatibility with older callers that still provide alpha.
alpha_rad = alpha_rad; %#ok<NASGU>

edgeSets = cell(nFN, 1);
for k = 1:nFN
    edgeSets{k} = zeros(0, 1, 'int32');
end

for k = 1:nFN
    ck = frameNodes(k, :);
    incidentElems = find(frameElems(:,1) == k | frameElems(:,2) == k);

    for e = incidentElems(:)'
        if frameElems(e,1) == k
            otherNode = frameElems(e,2);
        else
            otherNode = frameElems(e,1);
        end

        normal = frameNodes(otherNode, :) - ck;
        normalNorm = norm(normal);
        if normalNorm < eps('double')
            continue;
        end
        normal = normal / normalNorm;

        dp = boundaryCoords - ck;
        nearSection = vecnorm(dp, 2, 2) <= radialTol;
        if ~any(nearSection)
            continue;
        end

        planeDist = abs(dp * normal');
        snappedDist = round(planeDist(nearSection) / snapTol) * snapTol;
        uniqueDist = unique(snappedDist);
        minDist = uniqueDist(1);
        distSteps = diff(uniqueDist);
        distSteps = distSteps(distSteps > snapTol);
        if isempty(distSteps)
            planeTol = snapTol;
        else
            planeTol = 0.5 * min(distSteps);
        end

        onSection = nearSection & planeDist <= (minDist + planeTol + snapTol);
        edgeSets{k} = unique([edgeSets{k}; boundaryNodeIds(onSection)], 'stable');
    end
end
end

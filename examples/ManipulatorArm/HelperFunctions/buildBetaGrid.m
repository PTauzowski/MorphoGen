function B = buildBetaGrid(nJoints, deltaDeg, fixFirstJoint)
% buildBetaGrid  Full discrete grid of joint-angle combinations.
%
%   B = buildBetaGrid(nJoints, deltaDeg)
%
%   Returns an (N x nJoints) matrix where each row is one joint-angle
%   vector (in degrees) and N = (360/deltaDeg)^nJoints.
%
%   Each joint independently takes values {0, deltaDeg, 2*deltaDeg, ...
%   360-deltaDeg}.  The first joint is conventionally fixed at 0 (the
%   arm root has no free twist), so the effective grid size is
%   1 x (360/deltaDeg)^(nJoints-1) rows when fixFirstJoint=true (default).
%
%   Examples
%     B = buildBetaGrid(7, 90);   % 4^6 = 4096 rows  (first joint fixed)
%     B = buildBetaGrid(7, 45);   % 8^6 = 262144 rows
%
%   The returned matrix can be large for small deltaDeg.  Use
%   buildBetaGrid(nJoints, deltaDeg, false) to include the first joint.

    if nargin < 3
        fixFirstJoint = true;
    end

    vals = (0 : deltaDeg : 360 - deltaDeg)';   % column of angle values

    if fixFirstJoint
        % First joint fixed at 0; build grid over remaining joints.
        freeJoints = nJoints - 1;
        grids = cell(freeJoints, 1);
        for k = 1:freeJoints
            grids{k} = vals;
        end
        idx = buildFullGrid(grids);         % [N x freeJoints]
        B = [zeros(size(idx, 1), 1), idx];  % prepend zero first joint
    else
        grids = cell(nJoints, 1);
        for k = 1:nJoints
            grids{k} = vals;
        end
        B = buildFullGrid(grids);
    end
end

% -------------------------------------------------------------------------
function G = buildFullGrid(grids)
% Full Cartesian product of column vectors in the cell array grids.
% Returns [N x numel(grids)] matrix.
    n = numel(grids);
    sizes = cellfun(@numel, grids);
    N = prod(sizes);
    G = zeros(N, n);
    repEach = 1;
    repBlock = N;
    for k = 1:n
        repBlock = repBlock / sizes(k);
        col = repmat(repelem(grids{k}, repEach), repBlock, 1);
        G(:, k) = col;
        repEach = repEach * sizes(k);
    end
end

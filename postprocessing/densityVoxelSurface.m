function [V, F] = densityVoxelSurface(model, rho, threshold, nModules)
% densityVoxelSurface  Extract hard-threshold voxel surface from a hex density field.
%
%   [V, F] = densityVoxelSurface(model, rho)
%   [V, F] = densityVoxelSurface(model, rho, threshold, nModules)
%
%   Selects elements where rho > threshold, extracts exterior boundary quads,
%   and splits them into triangles. Result is a blocky voxelated surface.
%   Symmetric with densitySmoothSurface for use in plotTopology/exportTopology.

    if nargin < 3 || isempty(threshold), threshold = 0.5; end
    if nargin < 4 || isempty(nModules),  nModules  = 0;   end

    nodes = model.mesh.nodes;
    elems = model.mesh.elems;
    rho   = rho(:);

    if nModules > 0
        H        = model.halfSegmentNelems;
        nElemUse = min(nModules * 2 * H, size(elems, 1));
        elems    = elems(1:nElemUse, :);
        rho      = rho(1:nElemUse);
    end

    solidElems  = elems(rho > threshold, :);
    facePattern = model.fe.shapeFn.fcontours';
    quads       = buildBoundaryQuads(solidElems, facePattern);

    if isempty(quads)
        V = zeros(0, 3);
        F = zeros(0, 3);
        return;
    end

    tris           = [quads(:, [1 2 3]); quads(:, [1 3 4])];
    usedIds        = unique(tris(:));
    remap          = zeros(size(nodes, 1), 1);
    remap(usedIds) = 1:numel(usedIds);
    V = nodes(usedIds, :);
    F = remap(tris);
end

% -------------------------------------------------------------------------
function quads = buildBoundaryQuads(elems, facePattern)
    if isempty(elems)
        quads = zeros(0, size(facePattern, 2));
        return;
    end
    nElems   = size(elems, 1);
    allQuads = zeros(nElems * size(facePattern, 1), size(facePattern, 2), 'like', elems);
    for f = 1:size(facePattern, 1)
        rows = (f-1)*nElems + (1:nElems);
        allQuads(rows, :) = elems(:, facePattern(f, :));
    end
    sortedQuads = sort(allQuads, 2);
    [~, ~, ic]  = unique(sortedQuads, 'rows');
    counts      = accumarray(ic, 1);
    quads       = allQuads(counts(ic) == 1, :);
end

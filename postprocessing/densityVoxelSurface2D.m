function [V, F, faceLabels, boundary] = densityVoxelSurface2D(model, rho, threshold, nModules, labels)
% densityVoxelSurface2D  Hard-threshold 2D surface from a quad density field.
%
%   [V, F]                       — for patch rendering of solid quads
%   [V, F, faceLabels]           — one integer label per solid element face
%   [V, F, faceLabels, boundary] — boundary{i} is [Ki x 2] closed polygon
%
%   V  [Mx2]  unique node coordinates of solid elements
%   F  [Nx4]  solid element connectivity (compacted node indices)
%   boundary  cell array of closed polygon vertex arrays (for extrusion)

    if nargin < 3 || isempty(threshold), threshold = 0.5; end
    if nargin < 4 || isempty(nModules),  nModules  = 0;   end
    if nargin < 5,                        labels    = [];  end

    nodes  = model.mesh.nodes(:, 1:2);
    elems  = model.mesh.elems;
    rho    = rho(:);
    if ~isempty(labels), labels = labels(:); end

    if nModules > 0
        H        = model.halfSegmentNelems;
        nElemUse = min(nModules * 2 * H, size(elems, 1));
        elems    = elems(1:nElemUse, :);
        rho      = rho(1:nElemUse);
        if ~isempty(labels), labels = labels(1:nElemUse); end
    end

    solid   = rho > threshold;
    solidEl = elems(solid, :);

    if isempty(solidEl)
        V = zeros(0, 2);  F = zeros(0, size(elems,2));
        faceLabels = zeros(0, 1);
        boundary   = {};
        return;
    end

    % Compact node indices for patch rendering
    usedIds        = unique(solidEl(:));
    remap          = zeros(size(nodes, 1), 1);
    remap(usedIds) = 1:numel(usedIds);
    V = nodes(usedIds, :);
    F = remap(solidEl);

    % Per-face labels (one per solid element)
    if isempty(labels)
        faceLabels = zeros(nnz(solid), 1);
    else
        faceLabels = max(0, round(labels(solid)));
    end

    % Boundary polygons (for extrusion)
    if nargout > 3
        boundary = extractBoundaryPolygons2D(solidEl, nodes);
    end
end

% -------------------------------------------------------------------------
function boundary = extractBoundaryPolygons2D(solidElems, nodes)
    nEl = size(solidElems, 1);
    nV  = size(solidElems, 2);
    allEdges = zeros(nEl * nV, 2);
    for i = 1:nV
        j    = mod(i, nV) + 1;
        rows = (i-1)*nEl + (1:nEl);
        allEdges(rows, :) = solidElems(:, [i, j]);
    end
    sortedEdges = sort(allEdges, 2);
    [~, ~, ic]  = unique(sortedEdges, 'rows');
    counts      = accumarray(ic, 1);
    boundaryEdges = allEdges(counts(ic) == 1, :);
    boundary = chainEdges2D(boundaryEdges, nodes);
end

% -------------------------------------------------------------------------
function polygons = chainEdges2D(edges, nodes)
    polygons = {};
    if isempty(edges), return; end

    nEdges = size(edges, 1);
    used   = false(nEdges, 1);

    adj = cell(max(edges(:)), 1);
    for e = 1:nEdges
        adj{edges(e,1)} = [adj{edges(e,1)}, e];
        adj{edges(e,2)} = [adj{edges(e,2)}, e];
    end

    while any(~used)
        s        = find(~used, 1);
        used(s)  = true;
        chain    = [edges(s,1); edges(s,2)];

        while chain(end) ~= chain(1)
            cur    = chain(end);
            nbrs   = adj{cur};
            moved  = false;
            for e = nbrs
                if ~used(e)
                    used(e)    = true;
                    next       = edges(e,1);
                    if next == cur, next = edges(e,2); end
                    chain(end+1) = next; %#ok<AGROW>
                    moved = true;
                    break;
                end
            end
            if ~moved, break; end
        end

        polygons{end+1} = nodes(chain(1:end-1), :); %#ok<AGROW>
    end
end

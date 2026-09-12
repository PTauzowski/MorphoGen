function solid = removeDisconnectedComponents(solid, mesh, anchorElems)
% REMOVEDISCONNECTEDCOMPONENTS  Discard solid regions not connected to anchors.
%
%   solid = removeDisconnectedComponents(solid, mesh, anchorElems)
%
%   Flood-fills from anchorElems through the solid mask, keeping only the
%   reachable portion.  Two elements are considered connected when they
%   share at least one mesh node (node-adjacency), which is the most
%   generous connectivity and avoids splitting thin diagonal bridges.
%
%   Inputs:
%     solid        [nElems x 1] logical  binary topology mask
%     mesh         struct with .elems [nElems x nNodesPerElem]
%     anchorElems  [k x 1]  element indices to use as flood-fill seeds
%                           (typically from findAnchorElements)
%
%   Output:
%     solid  [nElems x 1] logical  mask with disconnected islands removed
%
%   If none of the anchorElems are solid, the mask is returned unchanged
%   with a warning (this prevents accidentally voiding everything when
%   anchor detection fails).

    solid = solid(:);
    elems    = mesh.elems;
    nElems   = size(elems, 1);
    nPerElem = size(elems, 2);
    nNodes   = max(elems(:));

    % Sparse node-to-element adjacency: N(node, elem) = 1
    N = sparse(elems(:), repmat((1:nElems)', nPerElem, 1), true, nNodes, nElems);

    seeds = anchorElems(solid(anchorElems));
    if isempty(seeds)
        warning('removeDisconnectedComponents: no solid anchor elements found; mask unchanged.');
        return;
    end

    visited = false(nElems, 1);
    visited(seeds) = true;
    frontier = seeds(:);

    while ~isempty(frontier)
        frontierNodes = find(any(N(:, frontier), 2));
        adjElems      = find(any(N(frontierNodes, :), 1))';
        newElems      = adjElems(~visited(adjElems) & solid(adjElems));
        visited(newElems) = true;
        frontier = newElems;
    end

    nRemoved = sum(solid) - sum(visited);
    if nRemoved > 0
        fprintf('removeDisconnectedComponents: removed %d disconnected solid elements (%.1f%% of solid).\n', ...
            nRemoved, 100 * nRemoved / sum(solid));
    end
    solid = visited;
end

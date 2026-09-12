function anchorElems = findAnchorElements(analyses, mesh)
% FINDANCHORELEMENTS  Element indices touching a support or load node in any config.
%
%   anchorElems = findAnchorElements(analyses, mesh)
%
%   Inputs:
%     analyses  {nConfigs x 1} cell of FEAnalysis objects, or a single object
%               Each must have .supports [nNodes x nDOFs] and
%               .Pnodal [nNodes x nDOFs] (inherited from FEAnalysis).
%     mesh      struct with .elems [nElems x nNodesPerElem] and
%               .nodes [nNodes x 3]  (full-arm mesh)
%
%   Output:
%     anchorElems  [k x 1]  element indices (into rows of mesh.elems)
%
%   Note: uses the full-arm mesh, so analyses must correspond to the
%   same mesh (same node ordering) that produced mesh.elems.

    if ~iscell(analyses)
        analyses = {analyses};
    end

    nNodes      = size(mesh.nodes, 1);
    anchorNodes = false(nNodes, 1);

    for k = 1:numel(analyses)
        a = analyses{k};
        anchorNodes = anchorNodes ...
            | any(a.supports ~= 0, 2) ...
            | any(a.Pnodal   ~= 0, 2);
    end

    anchorNodeIds = find(anchorNodes);
    anchorElems   = find(any(ismember(mesh.elems, anchorNodeIds), 2));
end

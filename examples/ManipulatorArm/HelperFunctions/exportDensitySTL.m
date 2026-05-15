function exportDensitySTL(model, rho, filename, threshold, nModules)
% exportDensitySTL  Write density-thresholded topology as binary STL.
%
%   exportDensitySTL(model, rho, filename)
%   exportDensitySTL(model, rho, filename, threshold)
%   exportDensitySTL(model, rho, filename, threshold, nModules)
%
%   Solid elements (rho > threshold) are extracted, their exterior boundary
%   quad faces are identified (faces shared by exactly one solid element),
%   split into triangles, and written as a binary STL file.
%
%   Inputs
%     model     - ManipulatorModel3D providing mesh.nodes, mesh.elems,
%                 fe.sf.fcontours, and halfSegmentNelems
%     rho       - [nElems x 1] element densities for the full arm
%     filename  - output file path (should end with .stl)
%     threshold - density threshold for solid/void classification (default 0.5)
%     nModules  - if given and > 0, export only the first nModules modules;
%                 each module spans 2*halfSegmentNelems elements (default: all)

    if nargin < 4 || isempty(threshold)
        threshold = 0.5;
    end
    if nargin < 5
        nModules = 0;   % 0 = export full arm
    end

    nodes = model.mesh.nodes;
    elems = model.mesh.elems;
    rho   = rho(:);

    % Optionally restrict to the first nModules modules.
    if nModules > 0
        H       = model.halfSegmentNelems;
        nElemUse = min(nModules * 2 * H, size(elems, 1));
        elems   = elems(1:nElemUse, :);
        rho     = rho(1:nElemUse);
    end

    selected = rho > threshold;
    if ~any(selected)
        warning('exportDensitySTL: no elements exceed threshold %.2f; %s not written.', ...
            threshold, filename);
        return;
    end
    solidElems = elems(selected, :);

    % Each hex8 face is a quad; fcontours is 4×6 (stored transposed) so
    % fcontours' is 6×4: each row lists the 4 local node indices of one face.
    facePattern = model.fe.sf.fcontours';   % [6 x 4]

    % Collect all quad faces and keep only boundary ones (appear exactly once).
    quads = buildBoundaryQuads(solidElems, facePattern);
    if isempty(quads)
        warning('exportDensitySTL: no boundary faces found; %s not written.', filename);
        return;
    end

    % Split each quad [a b c d] into triangles [a b c] and [a c d].
    tris = [quads(:, [1 2 3]); quads(:, [1 3 4])];

    % Compact node list to only nodes referenced by the surface triangles.
    usedIds        = unique(tris(:));
    remap          = zeros(size(nodes, 1), 1);
    remap(usedIds) = 1:numel(usedIds);
    V = nodes(usedIds, :);
    F = remap(tris);

    writeBinarySTL(filename, V, F);
end

% -------------------------------------------------------------------------
function quads = buildBoundaryQuads(elems, facePattern)
% Extract quad faces that are shared by exactly one solid element.
    nElems        = size(elems, 1);
    nFacesPerElem = size(facePattern, 1);
    nNodesPerFace = size(facePattern, 2);

    % Build the full [nElems*nFacesPerElem x nNodesPerFace] face table.
    allQuads = zeros(nElems * nFacesPerElem, nNodesPerFace, 'like', elems);
    for f = 1:nFacesPerElem
        localIdx = facePattern(f, :);                    % 1×4 local node ids
        rows     = (f-1)*nElems + 1 : f*nElems;
        allQuads(rows, :) = elems(:, localIdx);
    end

    % Boundary faces appear exactly once when sorted.
    sortedQuads  = sort(allQuads, 2);
    [~, ~, ic]   = unique(sortedQuads, 'rows');
    counts       = accumarray(ic, 1);
    quads        = allQuads(counts(ic) == 1, :);
end

% -------------------------------------------------------------------------
function writeBinarySTL(filename, V, F)
% Write binary STL from vertex array V [Mx3] and triangle array F [Nx3].
    nTri = size(F, 1);
    v1 = V(F(:,1), :);
    v2 = V(F(:,2), :);
    v3 = V(F(:,3), :);

    % Outward normals by right-hand rule on the face winding.
    normals = cross(v2 - v1, v3 - v1, 2);
    norms   = max(sqrt(sum(normals.^2, 2)), eps);
    normals = normals ./ norms;

    fid = fopen(filename, 'wb');
    if fid < 0
        error('exportDensitySTL: cannot open file for writing: %s', filename);
    end
    % 80-byte ASCII header.
    header = sprintf('%-80s', 'Binary STL exported by exportDensitySTL');
    fwrite(fid, header(1:80), 'char');
    fwrite(fid, uint32(nTri), 'uint32');
    for i = 1:nTri
        fwrite(fid, single(normals(i, :)), 'single');
        fwrite(fid, single(v1(i, :)),      'single');
        fwrite(fid, single(v2(i, :)),      'single');
        fwrite(fid, single(v3(i, :)),      'single');
        fwrite(fid, uint16(0),             'uint16');
    end
    fclose(fid);
end

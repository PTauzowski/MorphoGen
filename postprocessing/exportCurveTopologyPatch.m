function exportCurveTopologyPatch(model, rho, filename, threshold)
% exportCurveTopologyPatch  Save selected topology boundary faces as MAT data.
%
% The file contains vertices, faces, and selected element ids.  It is a light
% interchange format for later STL/mesh export without recomputing topology.

    if nargin < 4
        threshold = 0.5;
    end
    selected = rho(:) > threshold;
    elems = model.mesh.elems(selected, :);
    facePattern = model.fe.sf.fcontours';
    faces = elementFaces(elems, facePattern);
    faces = boundaryFacesOnly(faces);
    vertices = model.mesh.nodes; %#ok<NASGU>
    selectedElementIds = find(selected); %#ok<NASGU>
    save(filename, 'vertices', 'faces', 'selectedElementIds', 'threshold');
end

function faces = elementFaces(elems, facePattern)
    faces = zeros(size(elems, 1) * size(facePattern, 1), size(facePattern, 2));
    row = 1;
    for e = 1:size(elems, 1)
        ef = reshape(elems(e, facePattern(:)), size(facePattern));
        n = size(ef, 1);
        faces(row:row+n-1, :) = ef;
        row = row + n;
    end
end

function faces = boundaryFacesOnly(faces)
    if isempty(faces)
        return;
    end
    sortedFaces = sort(faces, 2);
    [~, ~, ic] = unique(sortedFaces, 'rows');
    counts = accumarray(ic, 1);
    faces = faces(counts(ic) == 1, :);
end

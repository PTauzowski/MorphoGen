function fig = plotCurveTopologyByOrigin(model, rho, originFull, originNames, titleText, threshold)
% plotCurveTopologyByOrigin  Real topology colored by dominant curve family.

    if nargin < 6
        threshold = 0.5;
    end

    fig = figure('Color', 'white', 'Visible', 'on');
    hold on; axis on; daspect([1 1 1]); view(45, 28);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(titleText, 'Interpreter', 'none');

    selected = rho(:) > threshold;
    colors = originColorMap();
    labels = ["mixed"; originNames(:)];
    handles = gobjects(0);
    legendLabels = strings(0, 1);

    for label = 0:numel(originNames)
        mask = selected & originFull(:) == label;
        if ~any(mask)
            continue;
        end
        color = colors(label + 1, :);
        h = plotSelectedBoundaryTopology(model, mask, color, 1.0);
        handles(end + 1, 1) = h; %#ok<AGROW>
        legendLabels(end + 1, 1) = labels(label + 1); %#ok<AGROW>
    end

    if ~isempty(handles)
        legend(handles, cellstr(legendLabels), 'Location', 'bestoutside', 'Interpreter', 'none');
    end
    camlight('headlight');
    lighting gouraud;
    axis equal;
end

function h = plotSelectedBoundaryTopology(model, selected, color, alphaValue)
    elems = model.mesh.elems(selected, :);
    if isempty(elems)
        h = gobjects(1);
        return;
    end
    facePattern = model.fe.sf.fcontours';
    faces = elementFaces(elems, facePattern);
    faces = boundaryFacesOnly(faces);
    h = patch('Vertices', model.mesh.nodes, 'Faces', faces, ...
        'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alphaValue, ...
        'BackFaceLighting', 'lit');
end

function colors = originColorMap()
    colors = [
        0.55 0.55 0.55;  % mixed
        0.88 0.10 0.10;  % helixPlus
        0.10 0.25 0.90;  % helixMinus
        0.10 0.62 0.24;  % axial
        0.95 0.55 0.05;  % bending
        0.55 0.20 0.75;  % ring
        0.95 0.85 0.05;  % jointRing
    ];
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

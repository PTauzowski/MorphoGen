function fig = plotCurveTopology(model, rho, titleText, mode)
% plotCurveTopology  Continuous or selected-element topology plot.

    if nargin < 4
        mode = "continuous";
    end
    mode = string(mode);

    fig = figure('Color', 'white', 'Visible', 'on');
    hold on; axis on; daspect([1 1 1]); view(45, 28);
    xlabel('x'); ylabel('y'); zlabel('z');

    switch mode
        case "continuous"
            plotElementDensityField(model, rho);
            caxis([0 1]);
            colormap(parula);
            cb = colorbar;
            cb.Label.String = 'density';
        case {"threshold", "binary"}
            plotTransparentShell(model, [0.78 0.78 0.78], 0.10);
            selected = rho(:) > 0.5;
            plotSelectedElements(model, selected, [0.95 0.76 0.12], 1.0);
        case {"real", "solid_topology", "void"}
            selected = rho(:) > 0.5;
            plotSelectedBoundaryTopology(model, selected, [0.95 0.76 0.12], 1.0);
        otherwise
            error('Unknown topology plot mode "%s".', mode);
    end

    title(titleText, 'Interpreter', 'none');
    axis equal;
end

function plotTransparentShell(model, color, alphaValue)
    elems = model.mesh.elems;
    facePattern = model.fe.sf.fcontours';
    faces = elementFaces(elems, facePattern);
    patch('Vertices', model.mesh.nodes, 'Faces', faces, ...
        'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alphaValue);
end

function plotSelectedElements(model, selected, color, alphaValue)
    elems = model.mesh.elems(selected, :);
    if isempty(elems)
        return;
    end
    facePattern = model.fe.sf.fcontours';
    faces = elementFaces(elems, facePattern);
    patch('Vertices', model.mesh.nodes, 'Faces', faces, ...
        'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alphaValue);
end

function plotSelectedBoundaryTopology(model, selected, color, alphaValue)
    elems = model.mesh.elems(selected, :);
    if isempty(elems)
        return;
    end
    facePattern = model.fe.sf.fcontours';
    faces = elementFaces(elems, facePattern);
    faces = boundaryFacesOnly(faces);
    patch('Vertices', model.mesh.nodes, 'Faces', faces, ...
        'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alphaValue, ...
        'BackFaceLighting', 'lit');
    camlight('headlight');
    lighting gouraud;
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

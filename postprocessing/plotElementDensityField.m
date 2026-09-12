function plotElementDensityField(model, x)
    hold on; axis on; daspect([1 1 1]); view(45, 35);
    colormap(parula); colorbar; caxis([0 1]);

    elems = model.mesh.elems;
    facePattern = model.fe.shapeFn.fcontours';
    faces = zeros(size(elems, 1) * size(facePattern, 1), size(facePattern, 2));
    faceColor = zeros(size(faces, 1), 1);

    row = 1;
    for e = 1:size(elems, 1)
        ef = reshape(elems(e, facePattern(:)), size(facePattern));
        n = size(ef, 1);
        faces(row:row+n-1, :) = ef;
        faceColor(row:row+n-1) = x(e);
        row = row + n;
    end

    patch('Vertices', model.mesh.nodes, 'Faces', faces, ...
        'FaceVertexCData', faceColor, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 1.0);
    xlabel('x'); ylabel('y'); zlabel('z');
end

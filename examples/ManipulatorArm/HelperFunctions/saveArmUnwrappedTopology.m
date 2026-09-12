function fig = saveArmUnwrappedTopology(model, x, resultRoot, stem, titleText, opts)
% saveArmUnwrappedTopology  Save unwrapped half-segment density/topology plots.
%
% Linked fields are shown as one reference half-segment.  Unlinked fields are
% shown as one panel per full-arm half-segment copy.

    if nargin < 6
        opts = struct();
    end
    if ~isfield(opts, 'threshold'), opts.threshold = []; end
    if ~isfield(opts, 'mode'), opts.mode = "auto"; end

    fig = plotArmUnwrappedTopology(model, x, titleText, opts);
    exportgraphics(fig, fullfile(resultRoot, [char(stem) '.png']), 'Resolution', 200);
    savefig(fig, fullfile(resultRoot, [char(stem) '.fig']));
    close(fig);
end

function fig = plotArmUnwrappedTopology(model, x, titleText, opts)
    x = x(:);
    H = model.halfSegmentNelems;
    nElems = size(model.mesh.elems, 1);
    assert(numel(x) == H || numel(x) == nElems, ...
        'Unwrapped plot density length must be H=%d or nElems=%d, got %d.', ...
        H, nElems, numel(x));

    mode = string(opts.mode);
    if mode == "auto"
        mode = inferMode(model, x);
    end

    [faces2d, verts2d] = referencePatchGeometry(model);

    fig = figure('Color', 'white', 'Visible', 'on', 'Name', titleText);
    if mode == "linked" || numel(x) == H
        ax = axes(fig);
        drawBlock(ax, faces2d, verts2d, localReferenceValues(model, x), opts.threshold);
        title(ax, sprintf('%s | linked reference', titleText), 'Interpreter', 'none');
        formatAxes(ax);
        addColorbarIfContinuous(ax, opts.threshold);
    else
        nCopies = nElems / H;
        nCols = min(4, nCopies);
        nRows = ceil(nCopies / nCols);
        tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
        for k = 1:nCopies
            ax = nexttile;
            ids = (k - 1) * H + (1:H);
            drawBlock(ax, faces2d, verts2d, x(ids), opts.threshold);
            title(ax, sprintf('copy %d', k), 'Interpreter', 'none');
            formatAxes(ax);
        end
        sgtitle(fig, sprintf('%s | unlinked full arm', titleText), 'Interpreter', 'none');
        if isempty(opts.threshold)
            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.Label.String = 'density';
        end
    end
end

function mode = inferMode(model, x)
    H = model.halfSegmentNelems;
    nElems = size(model.mesh.elems, 1);
    if numel(x) == H
        mode = "linked";
        return;
    end
    if mod(nElems, H) ~= 0
        mode = "allSegments";
        return;
    end

    ref = localReferenceValues(model, x);
    xLinked = model.segmentToArm(ref);
    if numel(xLinked) == numel(x) && max(abs(xLinked(:) - x(:))) < 1e-10
        mode = "linked";
    else
        mode = "allSegments";
    end
end

function ref = localReferenceValues(model, x)
    H = model.halfSegmentNelems;
    if numel(x) == H
        ref = x(:);
    else
        ref = x(1:H);
    end
end

function drawBlock(ax, faces2d, verts2d, values, threshold)
    values = values(:);
    if isempty(threshold)
        patch(ax, 'Faces', faces2d, 'Vertices', verts2d, ...
            'FaceVertexCData', repelem(values, 4), 'FaceColor', 'flat', ...
            'EdgeColor', 'none', 'FaceAlpha', 1.0);
        colormap(ax, parula);
        clim(ax, [0 1]);
    else
        selected = values > threshold;
        if any(selected)
            patch(ax, 'Faces', faces2d(selected, :), 'Vertices', verts2d, ...
                'FaceColor', [0.12 0.12 0.12], 'EdgeColor', 'none', ...
                'FaceAlpha', 1.0);
        end
    end
end

function formatAxes(ax)
    axis(ax, 'equal');
    axis(ax, 'tight');
    ax.Box = 'on';
    ax.TickDir = 'out';
    xlabel(ax, 'R theta [m]');
    ylabel(ax, 'local z [m]');
end

function addColorbarIfContinuous(ax, threshold)
    if isempty(threshold)
        cb = colorbar(ax);
        cb.Label.String = 'density';
    end
end

function [faces2d, verts2d] = referencePatchGeometry(model)
    H = model.halfSegmentNelems;
    elems = model.mesh.elems(1:H, :);
    nodes = model.mesh.nodes;
    R = model.R;

    verts2d = zeros(4 * H, 2);
    faces2d = reshape(1:(4 * H), 4, H).';
    for e = 1:H
        xyz = nodes(elems(e, :), :);
        theta = unwrap(atan2(xyz(:, 2), xyz(:, 1)));
        u = R * theta;
        z = xyz(:, 3);
        row = 4 * (e - 1) + 1;
        verts2d(row:row+3, :) = [
            min(u) min(z);
            max(u) min(z);
            max(u) max(z);
            min(u) max(z)
        ];
    end
    verts2d(:, 2) = verts2d(:, 2) - min(verts2d(:, 2));
end

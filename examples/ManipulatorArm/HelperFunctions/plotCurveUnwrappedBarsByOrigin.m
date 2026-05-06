function fig = plotCurveUnwrappedBarsByOrigin(model, rhoRef, originRef, originNames, params, titleText, threshold)
% plotCurveUnwrappedBarsByOrigin  Element-patch topology in unwrapped coordinates.

    if nargin < 7
        threshold = 0.5;
    end

    H = model.halfSegmentNelems;
    elems = model.mesh.elems(1:H, :);
    nodes = model.mesh.nodes;

    fig = figure('Color', 'white', 'Visible', 'on');
    hold on; axis equal tight;
    xlabel('R theta [m]');
    ylabel('local z [m]');
    title(titleText, 'Interpreter', 'none');

    colors = originColorMap();
    labels = ["mixed"; originNames(:)];
    handles = gobjects(0);
    legendLabels = strings(0, 1);

    for label = 0:numel(originNames)
        selectedElems = find(rhoRef(:) > threshold & originRef(:) == label);
        if isempty(selectedElems)
            continue;
        end
        [faces2d, verts2d] = unwrappedElementPatches(nodes, elems(selectedElems, :), model.R);
        h = patch('Faces', faces2d, 'Vertices', verts2d, ...
            'FaceColor', colors(label + 1, :), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
        handles(end + 1, 1) = h; %#ok<AGROW>
        legendLabels(end + 1, 1) = labels(label + 1); %#ok<AGROW>
    end

    if nargin >= 5 && ~isempty(params)
        xl = xlim;
        yl = ylim;
        overlayCurveFamiliesWithElementSize(model, params, xl(1), xl(2), yl(1), yl(2));
    end

    if ~isempty(handles)
        legend(handles, cellstr(legendLabels), 'Location', 'bestoutside', 'Interpreter', 'none');
    end
end

function [faces2d, verts2d] = unwrappedElementPatches(nodes, elems, R)
    nElem = size(elems, 1);
    verts2d = zeros(4 * nElem, 2);
    faces2d = reshape(1:(4 * nElem), 4, nElem).';

    for e = 1:nElem
        xyz = nodes(elems(e, :), :);
        theta = unwrap(atan2(xyz(:,2), xyz(:,1)));
        u = R * theta;
        z = xyz(:,3);
        uMin = min(u); uMax = max(u);
        zMin = min(z); zMax = max(z);
        row = 4 * (e - 1) + 1;
        verts2d(row:row+3, :) = [
            uMin zMin;
            uMax zMin;
            uMax zMax;
            uMin zMax
        ];
    end
end

function overlayCurveFamiliesWithElementSize(model, params, uMin, uMax, zMin, zMax)
    elemSize = estimateNominalElementSize(model);
    spacing = max(eps, params.spacingFactor * elemSize);
    angleRad = params.angleDeg * pi / 180;
    slope = tan(pi/2 - angleRad);
    phase = params.phaseFrac * spacing;
    zz = linspace(zMin, zMax, 200);
    kMin = floor((uMin - max(abs(slope * zz)) - phase) / spacing) - 1;
    kMax = ceil((uMax + max(abs(slope * zz)) - phase) / spacing) + 1;

    for k = kMin:kMax
        up = slope * zz + phase + k * spacing;
        um = -slope * zz - phase + k * spacing;
        plotClipped(up, zz, uMin, uMax, 'k-');
        plotClipped(um, zz, uMin, uMax, 'k-');
    end
end

function h = estimateNominalElementSize(model)
    h = (model.R - model.r) / max(1, model.resTh);
end

function plotClipped(u, z, uMin, uMax, style)
    mask = u >= uMin & u <= uMax;
    if nnz(mask) > 1
        plot(u(mask), z(mask), style, 'LineWidth', 0.6);
    end
end

function colors = originColorMap()
    colors = [
        0.55 0.55 0.55;
        0.88 0.10 0.10;
        0.10 0.25 0.90;
        0.10 0.62 0.24;
        0.95 0.55 0.05;
        0.55 0.20 0.75;
        0.95 0.85 0.05;
    ];
end

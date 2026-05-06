function fig = plotCurveUnwrappedTopology(model, rhoRef, params, titleText)
% plotCurveUnwrappedTopology  Plot reference half-segment density in u-z coordinates.

    H = model.halfSegmentNelems;
    elems = model.mesh.elems(1:H, :);
    nodes = model.mesh.nodes;
    centroids = squeeze(mean(reshape(nodes(elems', :), size(elems, 2), H, 3), 1));
    theta = atan2(centroids(:,2), centroids(:,1));
    u = model.R * theta;
    z = centroids(:,3);
    z = z - min(z);

    fig = figure('Color', 'white', 'Visible', 'on');
    scatter(u, z, 28, rhoRef(:), 'filled');
    hold on; axis tight;
    colormap(parula); caxis([0 1]);
    cb = colorbar; cb.Label.String = 'density';
    xlabel('R theta [m]');
    ylabel('local z [m]');
    title(titleText, 'Interpreter', 'none');

    if nargin >= 3 && ~isempty(params)
        overlayCurveFamilies(model, params, min(u), max(u), min(z), max(z));
    end
end


function overlayCurveFamilies(model, params, uMin, uMax, zMin, zMax)
    elemSize = 0.004;
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
        plotClipped(up, zz, uMin, uMax, 'r-');
        plotClipped(um, zz, uMin, uMax, 'b-');
    end
end

function plotClipped(u, z, uMin, uMax, style)
    mask = u >= uMin & u <= uMax;
    if nnz(mask) > 1
        plot(u(mask), z(mask), style, 'LineWidth', 0.8);
    end
end

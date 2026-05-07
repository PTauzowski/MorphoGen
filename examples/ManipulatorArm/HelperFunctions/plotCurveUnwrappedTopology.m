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
    legacyAngle   = curveFieldOrDefault(params, 'angleDeg', 45.0);
    anglePlusRad  = curveFieldOrDefault(params, 'anglePlusDeg',  legacyAngle) * pi / 180;
    angleMinusRad = curveFieldOrDefault(params, 'angleMinusDeg', legacyAngle) * pi / 180;
    slopePlus  = tan(pi/2 - anglePlusRad);
    slopeMinus = tan(pi/2 - angleMinusRad);
    phase = params.phaseFrac * spacing;
    zz = linspace(zMin, zMax, 200);
    maxDelta = max(abs(slopePlus), abs(slopeMinus)) * (zMax - zMin) + abs(phase);
    kMin = floor((uMin - uMax - maxDelta) / spacing) - 1;
    kMax = ceil((uMax - uMin + maxDelta) / spacing) + 1;

    for k = kMin:kMax
        up = slopePlus  * zz + phase + k * spacing;   % u - slopePlus*z  - phase = k*spacing
        um = -slopeMinus * zz - phase + k * spacing;  % u + slopeMinus*z + phase = k*spacing
        plotClipped(up, zz, uMin, uMax, 'r-');
        plotClipped(um, zz, uMin, uMax, 'b-');
    end
end

function v = curveFieldOrDefault(s, name, default)
    if isfield(s, name), v = s.(name); else, v = default; end
end

function plotClipped(u, z, uMin, uMax, style)
    mask = u >= uMin & u <= uMax;
    if nnz(mask) > 1
        plot(u(mask), z(mask), style, 'LineWidth', 0.8);
    end
end

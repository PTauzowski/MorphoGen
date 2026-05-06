function [rhoRef, rhoFull, info, fields] = buildCurveLinkedDensity(params, model, opts)
% buildCurveLinkedDensity  Generate linked full-arm density from curve params.
%
% The curve field is evaluated on the first half-segment element centroids.
% Helical, axial, bending-plane, periodic ring, and joint-ring families are
% combined. The result is rescaled to satisfy the requested full-arm volume
% fraction while keeping constant ring elements solid.

    H = model.halfSegmentNelems;
    elems = model.mesh.elems(1:H, :);
    nodes = model.mesh.nodes;
    centroids = squeeze(mean(reshape(nodes(elems', :), size(elems, 2), H, 3), 1));

    theta = atan2(centroids(:,2), centroids(:,1));
    z = centroids(:,3);
    z = z - min(z);
    zSpan = max(z) - min(z);
    if zSpan <= 0
        zSpan = model.R;
    end

    elemSize = localOpt(opts, 'elemSize', []);
    if isempty(elemSize)
        elemSize = medianElementSize(model, H);
    end

    spacing = max(eps, params.spacingFactor * elemSize);
    width = max(eps, params.widthFactor * elemSize);
    angleRad = params.angleDeg * pi / 180;
    slope = tan(pi/2 - angleRad); % du/dz in unwrapped u=R*theta coordinates
    phase = params.phaseFrac * spacing;

    u = model.R * theta;
    helixPlus = periodicLineDistance(u - slope * z - phase, spacing) ./ sqrt(1 + slope^2);
    helixMinus = periodicLineDistance(u + slope * z + phase, spacing) ./ sqrt(1 + slope^2);
    helixPlusField = paramOrDefault(params, 'helixPlusWeight', 1.0) * ...
        exp(-0.5 * (helixPlus / width).^2);
    helixMinusField = paramOrDefault(params, 'helixMinusWeight', 1.0) * ...
        exp(-0.5 * (helixMinus / width).^2);

    axialSpacing = max(spacing, 2 * elemSize);
    axialDist = periodicLineDistance(u - phase, axialSpacing);
    axialField = params.axialWeight * exp(-0.5 * (axialDist / width).^2);

    bendingDist = min(abs(wrapToPiLocal(theta)), abs(abs(wrapToPiLocal(theta)) - pi));
    bendingArcDist = model.R * bendingDist;
    bendingField = paramOrDefault(params, 'bendingWeight', 0.0) * ...
        exp(-0.5 * (bendingArcDist / width).^2);

    ringSpacingFactor = paramOrDefault(params, 'ringSpacingFactor', params.spacingFactor);
    ringSpacing = max(ringSpacingFactor * elemSize, 2 * elemSize);
    ringDist = periodicLineDistance(z - phase, ringSpacing);
    ringField = params.ringWeight * exp(-0.5 * (ringDist / width).^2);

    jointRingDist = min(z, zSpan - z);
    jointRingField = paramOrDefault(params, 'jointRingWeight', 0.0) * ...
        exp(-0.5 * (jointRingDist / width).^2);

    ridgeFields = [helixPlusField, helixMinusField, axialField, bendingField, ...
        ringField, jointRingField];
    raw = ridgeEnvelope(ridgeFields, 8);
    fieldsRef = struct();
    fieldsRef.helixPlus = helixPlusField;
    fieldsRef.helixMinus = helixMinusField;
    fieldsRef.axial = axialField;
    fieldsRef.bending = bendingField;
    fieldsRef.ring = ringField;
    fieldsRef.jointRing = jointRingField;
    fieldsRef.envelopeRaw = raw;

    raw = raw - min(raw);
    if max(raw) > 0
        raw = raw ./ max(raw);
    end
    fieldsRef.envelope = raw;

    rhoMin = localOpt(opts, 'rhoMin', 0.01);
    rhoRef = params.baseDensity + (1 - params.baseDensity) * raw;
    rhoRef = min(1, max(rhoMin, rhoRef));

    constRef = localOpt(opts, 'constRefElems', zeros(0, 1));
    rhoRef(constRef) = 1.0;

    rhoFull = model.segmentToArm(rhoRef);
    constFull = localOpt(opts, 'constFullElems', zeros(0, 1));
    rhoFull(constFull) = 1.0;

    targetVf = localOpt(opts, 'VolFrac', mean(rhoFull));
    [rhoFull, scaleInfo] = enforceDensityVolume(rhoFull, targetVf, rhoMin, constFull);

    % Pull the scaled first half-segment back for plotting/reporting.
    rhoRef = rhoFull(1:H);
    rhoRef(constRef) = 1.0;

    info = struct();
    info.linkedVolumeFraction = mean(rhoRef);
    info.fullVolumeFraction = mean(rhoFull);
    info.rawMin = min(raw);
    info.rawMax = max(raw);
    info.scale = scaleInfo.scale;
    info.targetVolumeFraction = targetVf;

    fields = struct();
    fields.reference = fieldsRef;
end

function d = periodicLineDistance(x, period)
    d = abs(mod(x + 0.5 * period, period) - 0.5 * period);
end

function a = wrapToPiLocal(a)
    a = mod(a + pi, 2*pi) - pi;
end

function value = paramOrDefault(params, name, defaultValue)
    if isfield(params, name)
        value = params.(name);
    else
        value = defaultValue;
    end
end

function raw = ridgeEnvelope(fields, q)
    fields = max(0, fields);
    if isempty(fields)
        raw = zeros(0, 1);
        return;
    end
    m = max(fields, [], 2);
    raw = zeros(size(m));
    active = m > 0;
    scaled = zeros(size(fields(active, :)));
    scaled(:, :) = fields(active, :) ./ m(active);
    raw(active) = m(active) .* sum(scaled .^ q, 2) .^ (1 / q);
end

function [rho, info] = enforceDensityVolume(rho, targetVf, rhoMin, fixedIds)
    rho = rho(:);
    fixed = false(size(rho));
    fixed(fixedIds) = true;
    free = ~fixed;
    rho(fixed) = 1.0;

    targetSum = targetVf * numel(rho);
    fixedSum = sum(rho(fixed));
    freeTarget = targetSum - fixedSum;
    if freeTarget <= rhoMin * nnz(free)
        rho(free) = rhoMin;
        scale = 0;
    elseif freeTarget >= nnz(free)
        rho(free) = 1.0;
        scale = Inf;
    else
        lo = rhoMin;
        base = max(0, rho(free) - lo);
        if sum(base) <= eps
            rho(free) = freeTarget / nnz(free);
            scale = 0;
        else
            scaleHi = max(1, (freeTarget - lo * nnz(free)) / sum(base));
            while sum(min(1, lo + scaleHi * base)) < freeTarget
                scaleHi = 2 * scaleHi;
            end
            scaleLo = 0;
            for it = 1:60
                scale = 0.5 * (scaleLo + scaleHi);
                trial = min(1, lo + scale * base);
                if sum(trial) < freeTarget
                    scaleLo = scale;
                else
                    scaleHi = scale;
                end
            end
            scale = scaleHi;
            rho(free) = min(1, lo + scale * base);
        end
    end
    info = struct('scale', scale);
end

function h = medianElementSize(model, H)
    elems = model.mesh.elems(1:H, :);
    nodes = model.mesh.nodes;
    centroids = squeeze(mean(reshape(nodes(elems', :), size(elems, 2), H, 3), 1));
    if size(centroids, 1) < 2
        h = model.R * 0.05;
        return;
    end
    sample = centroids(1:min(size(centroids, 1), 200), :);
    dmin = inf(size(sample, 1), 1);
    for i = 1:size(sample, 1)
        delta = sample - sample(i, :);
        d = sqrt(sum(delta.^2, 2));
        d(d <= 0) = inf;
        dmin(i) = min(d);
    end
    dmin = dmin(isfinite(dmin) & dmin > 0);
    h = median(dmin);
    if isempty(h) || ~isfinite(h) || h <= 0
        h = model.R * 0.05;
    end
end

function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name)
        value = s.(name);
    else
        value = defaultValue;
    end
end

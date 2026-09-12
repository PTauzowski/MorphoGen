function [polygons, polyLabels, info] = densitySmoothSurface2D(model, rho, threshold, nModules, opts)
% densitySmoothSurface2D  Smooth 2D isocurve from a quad density field.
%
%   [polygons, polyLabels, info] = densitySmoothSurface2D(model, rho)
%   [polygons, polyLabels, info] = densitySmoothSurface2D(model, rho, threshold, nModules, opts)
%
%   Extracts the rho=threshold boundary by:
%     1. Averaging element densities to mesh nodes
%     2. Optional Laplacian pre-smoothing of the nodal scalar field
%     3. Marching triangles on each quad element (split into 2 triangles)
%     4. Chaining boundary segments into closed polygons
%     5. Taubin polygon smoothing (shrink-free)
%
%   polygons   {K x 1} cell array of [Ni x 2] closed polygon vertex coords
%   polyLabels [K x 1] integer label per polygon (majority label of segments)
%   info       struct with smoothing parameters and polygon count
%
%   opts fields (all optional)
%     labels                [nElems x 1] integer label per element
%     nodalSmoothingIters   Laplacian pre-smooth iterations (default 1)
%     surfaceSmoothingIters Taubin polygon smooth iterations (default 12)
%     taubinLambda          positive Taubin step (default 0.5)
%     taubinMu              negative Taubin step (default -0.53)

    if nargin < 3 || isempty(threshold), threshold = 0.5;    end
    if nargin < 4 || isempty(nModules),  nModules  = 0;      end
    if nargin < 5 || isempty(opts),      opts      = struct(); end

    nodes  = model.mesh.nodes(:, 1:2);
    elems  = model.mesh.elems;
    rho    = rho(:);
    labels = localOpt(opts, 'labels', []);
    if ~isempty(labels), labels = labels(:); end

    if nModules > 0
        H        = model.halfSegmentNelems;
        nElemUse = min(nModules * 2 * H, size(elems, 1));
        elems    = elems(1:nElemUse, :);
        rho      = rho(1:nElemUse);
        if ~isempty(labels), labels = labels(1:nElemUse); end
    end

    if isempty(labels)
        labels = zeros(size(rho));
    else
        labels = max(0, round(labels(:)));
    end

    nodalSmoothingIters   = localOpt(opts, 'nodalSmoothingIters',   1);
    surfaceSmoothingIters = localOpt(opts, 'surfaceSmoothingIters', 12);
    taubinLambda          = localOpt(opts, 'taubinLambda',  0.50);
    taubinMu              = localOpt(opts, 'taubinMu',     -0.53);

    phi = elementToNodalDensity(elems, rho, size(nodes, 1));
    if nodalSmoothingIters > 0
        phi = smoothNodalScalar(elems, phi, nodalSmoothingIters);
    end

    [segments, segLabels] = marchingTriQuadMesh(nodes, elems, phi, threshold, labels);
    [polygons, polyLabels] = chainSegments(segments, segLabels);

    if surfaceSmoothingIters > 0
        for i = 1:numel(polygons)
            polygons{i} = taubinSmoothPolygon(polygons{i}, surfaceSmoothingIters, ...
                taubinLambda, taubinMu);
        end
    end

    info = struct( ...
        'nPolygons',             numel(polygons), ...
        'threshold',             threshold, ...
        'nodalSmoothingIters',   nodalSmoothingIters, ...
        'surfaceSmoothingIters', surfaceSmoothingIters);
end

% -------------------------------------------------------------------------
function phi = elementToNodalDensity(elems, rho, nNodes)
    nodeIds = elems(:);
    vals    = repmat(rho(:), size(elems, 2), 1);
    sumVals = accumarray(nodeIds, vals,  [nNodes, 1], @sum, 0);
    counts  = accumarray(nodeIds, 1,     [nNodes, 1], @sum, 0);
    phi     = zeros(nNodes, 1);
    used    = counts > 0;
    phi(used) = sumVals(used) ./ counts(used);
end

% -------------------------------------------------------------------------
function phi = smoothNodalScalar(elems, phi, nIter)
    nEl  = size(elems, 1);
    nV   = size(elems, 2);
    allEdges = zeros(nEl * nV, 2);
    for i = 1:nV
        j    = mod(i, nV) + 1;
        rows = (i-1)*nEl + (1:nEl);
        allEdges(rows, :) = elems(:, [i, j]);
    end
    edges = unique(sort(allEdges, 2), 'rows');
    n = numel(phi);
    for it = 1:nIter
        sumNbr = accumarray([edges(:,1); edges(:,2)], ...
            [phi(edges(:,2)); phi(edges(:,1))], [n, 1], @sum, 0);
        deg    = accumarray([edges(:,1); edges(:,2)], 1, [n, 1], @sum, 0);
        active = deg > 0;
        avg    = phi;
        avg(active) = sumNbr(active) ./ deg(active);
        phi(active) = 0.5 * phi(active) + 0.5 * avg(active);
    end
end

% -------------------------------------------------------------------------
function [segments, segLabels] = marchingTriQuadMesh(nodes, elems, phi, iso, labels)
    % Split each quad into 2 triangles and extract iso-boundary segments
    triPat    = [1 2 3; 1 3 4];
    segments  = zeros(0, 4);   % [x1 y1 x2 y2]
    segLabels = zeros(0, 1);

    for e = 1:size(elems, 1)
        eNodes = elems(e, :);
        eVals  = phi(eNodes);
        eX     = nodes(eNodes, :);
        if all(eVals >= iso) || all(eVals < iso), continue; end
        for ti = 1:2
            ids = triPat(ti, :);
            seg = marchingTri2D(eX(ids, :), eVals(ids), iso);
            if ~isempty(seg)
                segments(end+1, :)  = seg;        %#ok<AGROW>
                segLabels(end+1, 1) = labels(e);  %#ok<AGROW>
            end
        end
    end
end

% -------------------------------------------------------------------------
function seg = marchingTri2D(x, s, iso)
    % x [3x2], s [3x1] → seg [1x4] = [x1 y1 x2 y2] or []
    inside = s >= iso;
    nIn    = nnz(inside);
    seg    = [];
    if nIn == 0 || nIn == 3, return; end

    inIds  = find(inside);
    outIds = find(~inside);
    if nIn == 1
        a = inIds(1);   b = outIds(1);  c = outIds(2);
    else
        a = outIds(1);  b = inIds(1);   c = inIds(2);
    end
    p1  = interpIso2D(x(a,:), s(a), x(b,:), s(b), iso);
    p2  = interpIso2D(x(a,:), s(a), x(c,:), s(c), iso);
    seg = [p1, p2];
end

% -------------------------------------------------------------------------
function p = interpIso2D(x1, s1, x2, s2, iso)
    den = s2 - s1;
    if abs(den) <= eps, t = 0.5; else, t = (iso - s1) / den; end
    p = x1 + min(max(t, 0), 1) * (x2 - x1);
end

% -------------------------------------------------------------------------
function [polygons, polyLabels] = chainSegments(segments, segLabels)
    polygons   = {};
    polyLabels = zeros(0, 1);
    if isempty(segments), return; end

    % Merge near-duplicate endpoints using a coordinate key
    span = max(max(segments(:)) - min(segments(:)), 1);
    tol  = span * 1e-9;
    pts  = round([segments(:,1:2); segments(:,3:4)] ./ tol);
    [~, ia, ic] = unique(pts, 'rows');
    nSeg  = size(segments, 1);
    eA    = ic(1:nSeg);
    eB    = ic(nSeg+1:end);
    coords = [segments(:,1:2); segments(:,3:4)];
    coords = coords(ia, :);
    nPts  = size(coords, 1);

    adj = cell(nPts, 1);
    for s = 1:nSeg
        adj{eA(s)} = [adj{eA(s)}, s];
        adj{eB(s)} = [adj{eB(s)}, s];
    end

    used = false(nSeg, 1);
    while any(~used)
        s0       = find(~used, 1);
        used(s0) = true;
        chain    = [eA(s0); eB(s0)];
        labSum   = segLabels(s0);
        labCount = 1;

        while chain(end) ~= chain(1)
            cur   = chain(end);
            moved = false;
            for s = adj{cur}
                if ~used(s)
                    used(s)  = true;
                    labSum   = labSum + segLabels(s);
                    labCount = labCount + 1;
                    next     = eA(s);
                    if next == cur, next = eB(s); end
                    chain(end+1) = next; %#ok<AGROW>
                    moved = true;
                    break;
                end
            end
            if ~moved, break; end
        end

        polygons{end+1}   = coords(chain(1:end-1), :); %#ok<AGROW>
        polyLabels(end+1) = round(labSum / labCount);   %#ok<AGROW>
    end
    polyLabels = polyLabels(:);
end

% -------------------------------------------------------------------------
function poly = taubinSmoothPolygon(poly, nIter, lambda, mu)
    n = size(poly, 1);
    if n < 3, return; end
    prev = [n; (1:n-1)'];
    next = [(2:n)'; 1];
    for it = 1:nIter
        poly = polygonStep(poly, prev, next, lambda);
        poly = polygonStep(poly, prev, next, mu);
    end
end

function poly = polygonStep(poly, prev, next, factor)
    avg  = (poly(prev, :) + poly(next, :)) * 0.5;
    poly = poly + factor * (avg - poly);
end

% -------------------------------------------------------------------------
function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name), value = s.(name); else, value = defaultValue; end
end

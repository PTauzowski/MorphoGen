function info = exportDensitySmoothSTL(model, rho, filename, threshold, nModules, opts)
% exportDensitySmoothSTL  Write a smoothed density-isosurface STL.
%
%   exportDensitySmoothSTL(model, rho, filename)
%   exportDensitySmoothSTL(model, rho, filename, threshold, nModules, opts)
%
%   Unlike exportDensitySTL, this exporter does not threshold complete hex
%   elements. It first averages element densities to mesh nodes, extracts a
%   continuous rho=threshold surface with marching tetrahedra, caps portions
%   that reach the design-domain boundary, and optionally applies Taubin
%   surface smoothing.
%
%   This is intended for visualization/manufacturing-style geometry. Keep
%   exportDensitySTL for exact element-threshold diagnostics.

    if nargin < 4 || isempty(threshold)
        threshold = 0.5;
    end
    if nargin < 5 || isempty(nModules)
        nModules = 0;
    end
    if nargin < 6 || isempty(opts)
        opts = struct();
    end

    nodes = model.mesh.nodes;
    elems = model.mesh.elems;
    rho = rho(:);
    labels = localOpt(opts, 'labels', []);
    if ~isempty(labels)
        labels = labels(:);
    end

    if nModules > 0
        H = model.halfSegmentNelems;
        nElemUse = min(nModules * 2 * H, size(elems, 1));
        elems = elems(1:nElemUse, :);
        rho = rho(1:nElemUse);
        if ~isempty(labels)
            labels = labels(1:nElemUse);
        end
    end

    assert(numel(rho) == size(elems, 1), ...
        'Density length %d does not match element count %d.', numel(rho), size(elems, 1));
    if isempty(labels)
        labels = zeros(size(rho));
    else
        assert(numel(labels) == size(elems, 1), ...
            'Label length %d does not match element count %d.', numel(labels), size(elems, 1));
        labels = max(0, round(labels));
    end

    nodalSmoothingIters = localOpt(opts, 'nodalSmoothingIters', 1);
    surfaceSmoothingIters = localOpt(opts, 'surfaceSmoothingIters', 12);
    taubinLambda = localOpt(opts, 'taubinLambda', 0.45);
    taubinMu = localOpt(opts, 'taubinMu', -0.47);
    mergeTol = localOpt(opts, 'mergeTol', []);
    coloredObjPath = localOpt(opts, 'coloredObjPath', '');
    labelNames = string(localOpt(opts, 'labelNames', strings(0, 1)));

    phi = elementToNodalDensity(elems, rho, size(nodes, 1));
    if nodalSmoothingIters > 0
        phi = smoothNodalScalar(elems, phi, nodalSmoothingIters);
    end

    [V, F, faceLabels] = marchingTetHexMesh(nodes, elems, phi, threshold, labels);
    [Vb, Fb, faceLabelsB] = boundaryCaps(nodes, elems, phi, threshold, ...
        model.fe.sf.fcontours', labels);
    if ~isempty(Fb)
        Fb = Fb + size(V, 1);
        V = [V; Vb];
        F = [F; Fb];
        faceLabels = [faceLabels; faceLabelsB];
    end

    if isempty(F)
        warning('exportDensitySmoothSTL:noSurface', ...
            'No smooth isosurface found at threshold %.3g; %s not written.', ...
            threshold, filename);
        info = struct('nVertices', 0, 'nFaces', 0, 'threshold', threshold);
        return;
    end

    [V, F] = compactVertices(V, F, mergeTol);
    if surfaceSmoothingIters > 0
        V = taubinSmoothSurface(V, F, surfaceSmoothingIters, taubinLambda, taubinMu);
    end

    writeBinarySTL(filename, V, F);
    if ~isempty(coloredObjPath)
        writeColoredOBJ(coloredObjPath, V, F, faceLabels, labelNames, ...
            localOpt(opts, 'colors', []));
    end
    info = struct( ...
        'nVertices', size(V, 1), ...
        'nFaces', size(F, 1), ...
        'threshold', threshold, ...
        'nodalSmoothingIters', nodalSmoothingIters, ...
        'surfaceSmoothingIters', surfaceSmoothingIters);
end

% -------------------------------------------------------------------------
function phi = elementToNodalDensity(elems, rho, nNodes)
    nodeIds = elems(:);
    vals = repmat(rho(:), size(elems, 2), 1);
    sumVals = accumarray(nodeIds, vals, [nNodes, 1], @sum, 0);
    counts = accumarray(nodeIds, 1, [nNodes, 1], @sum, 0);
    phi = zeros(nNodes, 1);
    used = counts > 0;
    phi(used) = sumVals(used) ./ counts(used);
end

% -------------------------------------------------------------------------
function phi = smoothNodalScalar(elems, phi, nIter)
    edges = elementEdges(elems);
    n = numel(phi);
    for it = 1:nIter
        sumNbr = accumarray([edges(:,1); edges(:,2)], ...
            [phi(edges(:,2)); phi(edges(:,1))], [n, 1], @sum, 0);
        deg = accumarray([edges(:,1); edges(:,2)], 1, [n, 1], @sum, 0);
        active = deg > 0;
        avg = phi;
        avg(active) = sumNbr(active) ./ deg(active);
        phi(active) = 0.5 * phi(active) + 0.5 * avg(active);
    end
end

% -------------------------------------------------------------------------
function [V, F, faceLabels] = marchingTetHexMesh(nodes, elems, phi, iso, labels)
    tetPattern = [
        1 2 4 8
        1 4 3 8
        1 3 7 8
        1 7 5 8
        1 5 6 8
        1 6 2 8
    ];

    V = zeros(0, 3);
    F = zeros(0, 3);
    faceLabels = zeros(0, 1);
    for e = 1:size(elems, 1)
        eNodes = elems(e, :);
        eVals = phi(eNodes);
        if all(eVals >= iso) || all(eVals < iso)
            continue;
        end
        eX = nodes(eNodes, :);
        for ti = 1:size(tetPattern, 1)
            ids = tetPattern(ti, :);
            [vt, ft] = marchingTet(eX(ids, :), eVals(ids), iso);
            if isempty(ft)
                continue;
            end
            ft = ft + size(V, 1);
            V = [V; vt]; %#ok<AGROW>
            F = [F; ft]; %#ok<AGROW>
            faceLabels = [faceLabels; repmat(labels(e), size(ft, 1), 1)]; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
function [V, F] = marchingTet(x, s, iso)
    inside = s(:) >= iso;
    nIn = nnz(inside);
    V = zeros(0, 3);
    F = zeros(0, 3);
    if nIn == 0 || nIn == 4
        return;
    end

    inIds = find(inside);
    outIds = find(~inside);
    if nIn == 1 || nIn == 3
        if nIn == 1
            a = inIds(1);
            b = outIds(:)';
        else
            a = outIds(1);
            b = inIds(:)';
        end
        V = [
            interpIso(x(a,:), s(a), x(b(1),:), s(b(1)), iso)
            interpIso(x(a,:), s(a), x(b(2),:), s(b(2)), iso)
            interpIso(x(a,:), s(a), x(b(3),:), s(b(3)), iso)
        ];
        F = [1 2 3];
        return;
    end

    a = inIds(:)';
    b = outIds(:)';
    p11 = interpIso(x(a(1),:), s(a(1)), x(b(1),:), s(b(1)), iso);
    p12 = interpIso(x(a(1),:), s(a(1)), x(b(2),:), s(b(2)), iso);
    p21 = interpIso(x(a(2),:), s(a(2)), x(b(1),:), s(b(1)), iso);
    p22 = interpIso(x(a(2),:), s(a(2)), x(b(2),:), s(b(2)), iso);
    V = [p11; p12; p22; p21];
    F = [1 2 3; 1 3 4];
end

% -------------------------------------------------------------------------
function p = interpIso(x1, s1, x2, s2, iso)
    den = s2 - s1;
    if abs(den) <= eps
        t = 0.5;
    else
        t = (iso - s1) / den;
    end
    t = min(max(t, 0.0), 1.0);
    p = x1 + t * (x2 - x1);
end

% -------------------------------------------------------------------------
function [V, F, faceLabels] = boundaryCaps(nodes, elems, phi, iso, facePattern, labels)
    [quads, quadLabels] = buildBoundaryQuads(elems, facePattern, labels);
    V = zeros(0, 3);
    F = zeros(0, 3);
    faceLabels = zeros(0, 1);
    triPattern = [1 2 3; 1 3 4];
    for qi = 1:size(quads, 1)
        q = quads(qi, :);
        for ti = 1:2
            ids = q(triPattern(ti, :));
            poly = clipTriangleAboveIso(nodes(ids, :), phi(ids), iso);
            if size(poly, 1) < 3
                continue;
            end
            f0 = size(V, 1);
            V = [V; poly]; %#ok<AGROW>
            for k = 2:size(poly, 1)-1
                F = [F; f0 + [1 k k+1]]; %#ok<AGROW>
                faceLabels = [faceLabels; quadLabels(qi)]; %#ok<AGROW>
            end
        end
    end
end

% -------------------------------------------------------------------------
function poly = clipTriangleAboveIso(x, s, iso)
    poly = zeros(0, 3);
    for i = 1:3
        j = i + 1;
        if j > 3
            j = 1;
        end
        xi = x(i, :); xj = x(j, :);
        si = s(i);    sj = s(j);
        insideI = si >= iso;
        insideJ = sj >= iso;
        if insideI
            poly = [poly; xi]; %#ok<AGROW>
        end
        if insideI ~= insideJ
            poly = [poly; interpIso(xi, si, xj, sj, iso)]; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
function edges = elementEdges(elems)
    localEdges = [
        1 2; 2 4; 4 3; 3 1
        5 6; 6 8; 8 7; 7 5
        1 5; 2 6; 3 7; 4 8
    ];
    edges = zeros(size(elems, 1) * size(localEdges, 1), 2);
    for i = 1:size(localEdges, 1)
        rows = (i - 1) * size(elems, 1) + (1:size(elems, 1));
        edges(rows, :) = elems(:, localEdges(i, :));
    end
    edges = unique(sort(edges, 2), 'rows');
end

% -------------------------------------------------------------------------
function [quads, quadLabels] = buildBoundaryQuads(elems, facePattern, labels)
    nElems = size(elems, 1);
    allQuads = zeros(nElems * size(facePattern, 1), size(facePattern, 2), 'like', elems);
    allLabels = zeros(nElems * size(facePattern, 1), 1);
    for f = 1:size(facePattern, 1)
        rows = (f - 1) * nElems + (1:nElems);
        allQuads(rows, :) = elems(:, facePattern(f, :));
        allLabels(rows) = labels(:);
    end
    sortedQuads = sort(allQuads, 2);
    [~, ~, ic] = unique(sortedQuads, 'rows');
    counts = accumarray(ic, 1);
    keep = counts(ic) == 1;
    quads = allQuads(keep, :);
    quadLabels = allLabels(keep);
end

% -------------------------------------------------------------------------
function [V, F] = compactVertices(V, F, tol)
    if isempty(tol)
        span = max(max(V, [], 1) - min(V, [], 1));
        tol = max(span, 1) * 1.0e-10;
    end
    key = round(V ./ tol);
    [~, ia, ic] = unique(key, 'rows');
    V = V(ia, :);
    F = ic(F);
end

% -------------------------------------------------------------------------
function V = taubinSmoothSurface(V, F, nIter, lambda, mu)
    edges = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    for it = 1:nIter
        V = smoothStep(V, edges, lambda);
        V = smoothStep(V, edges, mu);
    end
end

function V = smoothStep(V, edges, factor)
    n = size(V, 1);
    sumNbr = zeros(n, 3);
    for d = 1:3
        sumNbr(:, d) = accumarray([edges(:,1); edges(:,2)], ...
            [V(edges(:,2), d); V(edges(:,1), d)], [n, 1], @sum, 0);
    end
    deg = accumarray([edges(:,1); edges(:,2)], 1, [n, 1], @sum, 0);
    active = deg > 0;
    avg = V;
    avg(active, :) = sumNbr(active, :) ./ deg(active);
    V(active, :) = V(active, :) + factor * (avg(active, :) - V(active, :));
end

% -------------------------------------------------------------------------
function writeBinarySTL(filename, V, F)
    nTri = size(F, 1);
    v1 = V(F(:,1), :);
    v2 = V(F(:,2), :);
    v3 = V(F(:,3), :);
    normals = cross(v2 - v1, v3 - v1, 2);
    norms = max(sqrt(sum(normals.^2, 2)), eps);
    normals = normals ./ norms;

    fid = fopen(filename, 'wb');
    if fid < 0
        error('exportDensitySmoothSTL:cannotOpenFile', ...
            'Cannot open file for writing: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));
    header = sprintf('%-80s', 'Binary STL exported by exportDensitySmoothSTL');
    fwrite(fid, header(1:80), 'char');
    fwrite(fid, uint32(nTri), 'uint32');
    for i = 1:nTri
        fwrite(fid, single(normals(i, :)), 'single');
        fwrite(fid, single(v1(i, :)), 'single');
        fwrite(fid, single(v2(i, :)), 'single');
        fwrite(fid, single(v3(i, :)), 'single');
        fwrite(fid, uint16(0), 'uint16');
    end
    clear cleaner;
end

% -------------------------------------------------------------------------
function writeColoredOBJ(filename, V, F, faceLabels, labelNames, colors)
    faceLabels = max(0, round(faceLabels(:)));
    maxLabel = max(faceLabels);
    names = ["mixed"; labelNames(:)];
    if maxLabel + 1 > numel(names)
        extra = "label" + string(numel(names):maxLabel);
        names = [names; extra(:)];
    end
    if isempty(colors)
        colors = defaultOriginColors(numel(names));
    elseif size(colors, 1) < numel(names)
        colors = [colors; lines(numel(names) - size(colors, 1))];
    end

    [objPath, mtlName] = normalizeObjPath(filename);
    mtlPath = fullfile(fileparts(objPath), mtlName);
    writeMtl(mtlPath, names, colors);

    fid = fopen(objPath, 'w');
    if fid < 0
        error('exportDensitySmoothSTL:cannotOpenObj', ...
            'Cannot open OBJ file for writing: %s', objPath);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# Smooth colored OBJ exported by exportDensitySmoothSTL\n');
    fprintf(fid, 'mtllib %s\n', mtlName);
    for i = 1:size(V, 1)
        fprintf(fid, 'v %.9g %.9g %.9g\n', V(i, 1), V(i, 2), V(i, 3));
    end
    fprintf(fid, '\n');

    labelsUsed = unique(faceLabels);
    for i = 1:numel(labelsUsed)
        lab = labelsUsed(i);
        matName = materialName(names(lab + 1), lab);
        fprintf(fid, 'g %s\n', matName);
        fprintf(fid, 'usemtl %s\n', matName);
        idx = find(faceLabels == lab);
        for k = 1:numel(idx)
            tri = F(idx(k), :);
            fprintf(fid, 'f %d %d %d\n', tri(1), tri(2), tri(3));
        end
        fprintf(fid, '\n');
    end
end

% -------------------------------------------------------------------------
function [objPath, mtlName] = normalizeObjPath(filename)
    [folder, stem, ext] = fileparts(filename);
    if isempty(ext)
        ext = '.obj';
    end
    objPath = fullfile(folder, stem + string(ext));
    mtlName = stem + ".mtl";
end

% -------------------------------------------------------------------------
function writeMtl(path, names, colors)
    fid = fopen(path, 'w');
    if fid < 0
        error('exportDensitySmoothSTL:cannotOpenMtl', ...
            'Cannot open MTL file for writing: %s', path);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# Materials exported by exportDensitySmoothSTL\n');
    for i = 1:numel(names)
        matName = materialName(names(i), i - 1);
        c = colors(i, :);
        fprintf(fid, 'newmtl %s\n', matName);
        fprintf(fid, 'Kd %.6f %.6f %.6f\n', c(1), c(2), c(3));
        fprintf(fid, 'Ka %.6f %.6f %.6f\n', 0.25*c(1), 0.25*c(2), 0.25*c(3));
        fprintf(fid, 'Ks 0.100000 0.100000 0.100000\n');
        fprintf(fid, 'Ns 24.000000\n');
        fprintf(fid, 'd 1.000000\n\n');
    end
end

% -------------------------------------------------------------------------
function matName = materialName(name, label)
    raw = char(name);
    raw = regexprep(raw, '[^A-Za-z0-9_]+', '_');
    if isempty(raw)
        raw = sprintf('label%d', label);
    end
    matName = sprintf('mat_%02d_%s', label, raw);
end

% -------------------------------------------------------------------------
function colors = defaultOriginColors(n)
    base = [
        0.55 0.55 0.55
        0.88 0.10 0.10
        0.10 0.25 0.90
        0.10 0.62 0.24
        0.95 0.55 0.05
        0.55 0.20 0.75
        0.95 0.85 0.05
        0.00 0.70 0.85
    ];
    if n <= size(base, 1)
        colors = base(1:n, :);
    else
        colors = [base; lines(n - size(base, 1))];
    end
end

% -------------------------------------------------------------------------
function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name)
        value = s.(name);
    else
        value = defaultValue;
    end
end

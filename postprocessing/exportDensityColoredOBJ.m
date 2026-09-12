function info = exportDensityColoredOBJ(model, rho, labels, labelNames, filename, threshold, nModules, opts)
% exportDensityColoredOBJ  Export thresholded topology as colored OBJ groups.
%
%   exportDensityColoredOBJ(model, rho, labels, labelNames, filename)
%
%   labels is one integer label per element. Label 0 is treated as "mixed";
%   labels 1..N correspond to labelNames. The exporter writes filename.obj
%   plus a sibling filename.mtl and assigns one material per label.
%
%   This is for visualization. STL does not have portable color/material
%   support, so OBJ+MTL is used instead.

    if nargin < 6 || isempty(threshold)
        threshold = 0.5;
    end
    if nargin < 7 || isempty(nModules)
        nModules = 0;
    end
    if nargin < 8 || isempty(opts)
        opts = struct();
    end

    nodes = model.mesh.nodes;
    elems = model.mesh.elems;
    rho = rho(:);
    labels = labels(:);

    if nModules > 0
        H = model.halfSegmentNelems;
        nElemUse = min(nModules * 2 * H, size(elems, 1));
        elems = elems(1:nElemUse, :);
        rho = rho(1:nElemUse);
        labels = labels(1:nElemUse);
    end

    assert(numel(rho) == size(elems, 1), ...
        'Density length %d does not match element count %d.', numel(rho), size(elems, 1));
    assert(numel(labels) == size(elems, 1), ...
        'Label length %d does not match element count %d.', numel(labels), size(elems, 1));

    labels = max(0, round(labels));
    labelNames = string(labelNames(:));
    names = ["mixed"; labelNames];
    maxLabel = max(labels);
    if maxLabel + 1 > numel(names)
        extra = "label" + string(numel(names):maxLabel);
        names = [names; extra(:)];
    end

    colors = localOpt(opts, 'colors', defaultOriginColors(numel(names)));
    if size(colors, 1) < numel(names)
        colors = [colors; lines(numel(names) - size(colors, 1))];
    end

    selected = rho > threshold;
    if ~any(selected)
        warning('exportDensityColoredOBJ:noElements', ...
            'No elements exceed threshold %.3g; %s not written.', threshold, filename);
        info = struct('nVertices', 0, 'nFaces', 0, 'nGroups', 0);
        return;
    end

    labelsUsed = unique(labels(selected));
    labelFaces = cell(numel(labelsUsed), 1);
    labelFaceLabels = zeros(numel(labelsUsed), 1);
    totalFaces = 0;

    facePattern = model.fe.shapeFn.fcontours';
    for i = 1:numel(labelsUsed)
        lab = labelsUsed(i);
        mask = selected & labels == lab;
        quads = buildBoundaryQuads(elems(mask, :), facePattern);
        tris = [quads(:, [1 2 3]); quads(:, [1 3 4])];
        labelFaces{i} = tris;
        labelFaceLabels(i) = lab;
        totalFaces = totalFaces + size(tris, 1);
    end

    usedIds = unique(vertcat(labelFaces{:}));
    remap = zeros(size(nodes, 1), 1);
    remap(usedIds) = 1:numel(usedIds);
    V = nodes(usedIds, :);

    [objPath, mtlName] = normalizeObjPath(filename);
    mtlPath = fullfile(fileparts(objPath), mtlName);
    writeMtl(mtlPath, names, colors);
    writeObj(objPath, mtlName, V, labelFaces, labelFaceLabels, remap, names);

    info = struct( ...
        'nVertices', size(V, 1), ...
        'nFaces', totalFaces, ...
        'nGroups', numel(labelsUsed), ...
        'objPath', string(objPath), ...
        'mtlPath', string(mtlPath));
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
        error('exportDensityColoredOBJ:cannotOpenMtl', ...
            'Cannot open MTL file for writing: %s', path);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# Materials exported by exportDensityColoredOBJ\n');
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
function writeObj(path, mtlName, V, labelFaces, labelFaceLabels, remap, names)
    fid = fopen(path, 'w');
    if fid < 0
        error('exportDensityColoredOBJ:cannotOpenObj', ...
            'Cannot open OBJ file for writing: %s', path);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# OBJ exported by exportDensityColoredOBJ\n');
    fprintf(fid, 'mtllib %s\n', mtlName);
    for i = 1:size(V, 1)
        fprintf(fid, 'v %.9g %.9g %.9g\n', V(i, 1), V(i, 2), V(i, 3));
    end
    fprintf(fid, '\n');

    for i = 1:numel(labelFaces)
        lab = labelFaceLabels(i);
        matName = materialName(names(lab + 1), lab);
        fprintf(fid, 'g %s\n', matName);
        fprintf(fid, 'usemtl %s\n', matName);
        tris = remap(labelFaces{i});
        for f = 1:size(tris, 1)
            fprintf(fid, 'f %d %d %d\n', tris(f, 1), tris(f, 2), tris(f, 3));
        end
        fprintf(fid, '\n');
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
function quads = buildBoundaryQuads(elems, facePattern)
    if isempty(elems)
        quads = zeros(0, size(facePattern, 2));
        return;
    end
    nElems = size(elems, 1);
    allQuads = zeros(nElems * size(facePattern, 1), size(facePattern, 2), 'like', elems);
    for f = 1:size(facePattern, 1)
        rows = (f - 1) * nElems + (1:nElems);
        allQuads(rows, :) = elems(:, facePattern(f, :));
    end
    sortedQuads = sort(allQuads, 2);
    [~, ~, ic] = unique(sortedQuads, 'rows');
    counts = accumarray(ic, 1);
    quads = allQuads(counts(ic) == 1, :);
end

% -------------------------------------------------------------------------
function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name)
        value = s.(name);
    else
        value = defaultValue;
    end
end

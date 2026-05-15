function info = exportDensitySmoothSTL(model, rho, filename, threshold, nModules, opts)
% exportDensitySmoothSTL  Write a smoothed density-isosurface STL.
%
%   exportDensitySmoothSTL(model, rho, filename)
%   exportDensitySmoothSTL(model, rho, filename, threshold, nModules, opts)
%
%   Unlike exportDensitySTL, this exporter extracts a continuous isosurface
%   via marching tetrahedra and applies Taubin surface smoothing.
%   See densitySmoothSurface for the full list of smoothing opts.
%
%   Additional opts
%     coloredObjPath   path to write a colored OBJ/MTL alongside the STL
%     labelNames       string array of label names (used with opts.labels)
%     colors           [nLabels x 3] RGB colors

    if nargin < 4 || isempty(threshold), threshold = 0.5;    end
    if nargin < 5 || isempty(nModules),  nModules  = 0;      end
    if nargin < 6 || isempty(opts),      opts      = struct(); end

    coloredObjPath = localOpt(opts, 'coloredObjPath', '');
    labelNames     = string(localOpt(opts, 'labelNames', strings(0, 1)));
    colors         = localOpt(opts, 'colors', []);

    [V, F, faceLabels, geoInfo] = densitySmoothSurface(model, rho, threshold, nModules, opts);

    if isempty(F)
        warning('exportDensitySmoothSTL:noSurface', ...
            'No smooth isosurface found at threshold %.3g; %s not written.', ...
            threshold, filename);
        info = struct('nVertices', 0, 'nFaces', 0, 'threshold', threshold);
        return;
    end

    writeBinarySTL(filename, V, F);

    if ~isempty(coloredObjPath)
        writeColoredOBJ(coloredObjPath, V, F, faceLabels, labelNames, colors);
    end

    info = geoInfo;
end

% -------------------------------------------------------------------------
function writeBinarySTL(filename, V, F)
    nTri = size(F, 1);
    v1 = V(F(:,1), :);  v2 = V(F(:,2), :);  v3 = V(F(:,3), :);
    normals = cross(v2 - v1, v3 - v1, 2);
    normals = normals ./ max(sqrt(sum(normals.^2, 2)), eps);

    fid = fopen(filename, 'wb');
    if fid < 0
        error('exportDensitySmoothSTL:cannotOpenFile', ...
            'Cannot open file for writing: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, sprintf('%-80s', 'Binary STL exported by exportDensitySmoothSTL'), 'char');
    fwrite(fid, uint32(nTri), 'uint32');
    for i = 1:nTri
        fwrite(fid, single(normals(i,:)), 'single');
        fwrite(fid, single(v1(i,:)),      'single');
        fwrite(fid, single(v2(i,:)),      'single');
        fwrite(fid, single(v3(i,:)),      'single');
        fwrite(fid, uint16(0),            'uint16');
    end
end

% -------------------------------------------------------------------------
function writeColoredOBJ(filename, V, F, faceLabels, labelNames, colors)
    faceLabels = max(0, round(faceLabels(:)));
    names      = ["mixed"; labelNames(:)];
    maxLabel   = max(faceLabels);
    if maxLabel + 1 > numel(names)
        names = [names; "label" + string(numel(names):maxLabel)'];
    end
    if isempty(colors)
        colors = defaultOriginColors(numel(names));
    elseif size(colors, 1) < numel(names)
        colors = [colors; lines(numel(names) - size(colors, 1))];
    end

    [objPath, mtlName] = normalizeObjPath(filename);
    writeMtl(fullfile(fileparts(objPath), mtlName), names, colors);

    fid = fopen(objPath, 'w');
    if fid < 0
        error('exportDensitySmoothSTL:cannotOpenObj', ...
            'Cannot open OBJ file for writing: %s', objPath);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# Smooth colored OBJ exported by exportDensitySmoothSTL\n');
    fprintf(fid, 'mtllib %s\n', mtlName);
    for i = 1:size(V, 1)
        fprintf(fid, 'v %.9g %.9g %.9g\n', V(i,1), V(i,2), V(i,3));
    end
    fprintf(fid, '\n');
    for lab = unique(faceLabels)'
        matName = materialName(names(lab + 1), lab);
        fprintf(fid, 'g %s\nusemtl %s\n', matName, matName);
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
    if isempty(ext), ext = '.obj'; end
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
        fprintf(fid, 'newmtl %s\nKd %.6f %.6f %.6f\nKa %.6f %.6f %.6f\n', ...
            matName, c(1), c(2), c(3), 0.25*c(1), 0.25*c(2), 0.25*c(3));
        fprintf(fid, 'Ks 0.1 0.1 0.1\nNs 24\nd 1\n\n');
    end
end

% -------------------------------------------------------------------------
function matName = materialName(name, label)
    raw = regexprep(char(name), '[^A-Za-z0-9_]+', '_');
    if isempty(raw), raw = sprintf('label%d', label); end
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
    if isstruct(s) && isfield(s, name), value = s.(name); else, value = defaultValue; end
end

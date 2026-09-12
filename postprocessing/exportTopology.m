function info = exportTopology(model, rho, resultRoot, opts)
% exportTopology  Write topology geometry to external file formats (2D or 3D).
%
%   exportTopology(model, rho, resultRoot)
%   exportTopology(model, rho, resultRoot, opts)
%   info = exportTopology(...)   % returns struct of written file paths
%
%   Auto-detects 2D (quad mesh) vs 3D (hex mesh) from model.mesh.nodes.
%   For 2D, STL and OBJ are produced by extruding the boundary profile in Z.
%
%   Supported formats
%     "stl"  Binary STL
%     "obj"  Colored OBJ + MTL
%     "mat"  MATLAB patch data for later reimport
%
%   Format opts
%     formats          string | string array, subset of ["stl","obj","mat"]
%                      (default "stl")
%     smoothed         true | false (default true)
%     colored          true | false — per-family coloring for OBJ (default false)
%     threshold        iso-value (default 0.5)
%     nModules         restrict to first N modules (default: all)
%     filenameStem     base name for all written files (default 'topology')
%
%   2D-specific opts
%     extrusionThickness  Z-depth for extruded STL/OBJ (default 1.0)
%
%   Coloring opts
%     labels           [nElems x 1] integer label per element
%     labelNames       string array of label names
%     colors           [nLabels x 3] RGB
%
%   Smoothing opts (passed through to geometry extraction)
%     nodalSmoothingIters    (default 1)
%     surfaceSmoothingIters  (default 12)
%     taubinLambda / taubinMu

    if nargin < 4 || isempty(opts), opts = struct(); end

    formats      = string(localOpt(opts, 'formats',      "stl"));
    smoothed     = localOpt(opts, 'smoothed',     true);
    colored      = localOpt(opts, 'colored',      false) | ~isempty(localOpt(opts, 'labels', []));
    threshold    = localOpt(opts, 'threshold',    0.5);
    nModules     = localOpt(opts, 'nModules',     0);
    filenameStem = string(localOpt(opts, 'filenameStem', 'topology'));
    labelNames   = string(localOpt(opts, 'labelNames',   strings(0, 1)));
    colors       = localOpt(opts, 'colors',   []);
    labels       = localOpt(opts, 'labels',   []);

    wantSTL = any(formats == "stl");
    wantOBJ = any(formats == "obj");
    wantMAT = any(formats == "mat");

    if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

    stlPath = fullfile(resultRoot, filenameStem + ".stl");
    objPath = fullfile(resultRoot, filenameStem + "_colored.obj");
    matPath = fullfile(resultRoot, filenameStem + ".mat");

    info = struct('stl', '', 'obj', '', 'mtl', '', 'mat', '');

    if is2D(model)
        info = exportTopology2D(model, rho, resultRoot, opts, ...
            wantSTL, wantOBJ, wantMAT, smoothed, colored, threshold, nModules, ...
            filenameStem, stlPath, objPath, matPath, labelNames, colors, labels);
    else
        info = exportTopology3D(model, rho, resultRoot, opts, ...
            wantSTL, wantOBJ, wantMAT, smoothed, colored, threshold, nModules, ...
            filenameStem, stlPath, objPath, matPath, labelNames, colors, labels, info);
    end
end

% =========================================================================
function info = exportTopology2D(model, rho, ~, opts, ...
        wantSTL, wantOBJ, wantMAT, smoothed, colored, threshold, nModules, ...
        filenameStem, stlPath, objPath, matPath, labelNames, colors, labels)

    thickness = localOpt(opts, 'extrusionThickness', 1.0);
    info      = struct('stl', '', 'obj', '', 'mtl', '', 'mat', '');

    if wantSTL || wantOBJ
        % Extract 2D boundary polygons
        if smoothed
            smoothOpts           = opts;
            smoothOpts.labels    = labels;
            [polygons, polyLabels] = densitySmoothSurface2D( ...
                model, rho, threshold, nModules, smoothOpts);
        else
            [~, ~, faceLabels2D, boundary] = densityVoxelSurface2D( ...
                model, rho, threshold, nModules, labels);
            polygons   = boundary;
            polyLabels = zeros(numel(polygons), 1);
            if ~isempty(faceLabels2D) && ~isempty(polygons)
                polyLabels(:) = mode(double(faceLabels2D));
            end
        end

        if isempty(polygons)
            warning('exportTopology:noGeometry2D', ...
                'No solid region found at threshold %.3g.', threshold);
        else
            % Extrude to 3D
            [V3, F3, fl3] = extrudeProfile2D(polygons, polyLabels, thickness);

            if wantSTL
                writeBinarySTL(stlPath, V3, F3);
                info.stl = stlPath;
            end

            if wantOBJ && colored
                writeColoredOBJ(objPath, V3, F3, fl3, labelNames, colors, filenameStem);
                info.obj = objPath;
                info.mtl = strrep(objPath, '.obj', '.mtl');
            elseif wantOBJ
                writeColoredOBJ(objPath, V3, F3, zeros(size(F3,1),1), ...
                    strings(0,1), [], filenameStem);
                info.obj = objPath;
                info.mtl = strrep(objPath, '.obj', '.mtl');
            end
        end
    end

    if wantMAT
        exportCurveTopologyPatch(model, rho, matPath, threshold);
        info.mat = matPath;
    end
end

% =========================================================================
function info = exportTopology3D(model, rho, ~, opts, ...
        wantSTL, wantOBJ, wantMAT, smoothed, colored, threshold, nModules, ...
        filenameStem, stlPath, objPath, matPath, labelNames, colors, labels, info)

    % --- STL and smooth colored OBJ (share one geometry pass) ------------
    if wantSTL || (wantOBJ && smoothed)
        exportOpts            = opts;
        exportOpts.labels     = labels;
        exportOpts.labelNames = labelNames;
        exportOpts.colors     = colors;

        if wantOBJ && smoothed && colored
            exportOpts.coloredObjPath = objPath;
        end

        actualStlPath = stlPath;
        deleteTempSTL = false;
        if ~wantSTL
            actualStlPath = [tempname '.stl'];
            deleteTempSTL = true;
        end

        if smoothed
            exportDensitySmoothSTL(model, rho, actualStlPath, threshold, nModules, exportOpts);
        else
            exportDensitySTL(model, rho, actualStlPath, threshold, nModules);
        end

        if wantSTL
            info.stl = actualStlPath;
        elseif deleteTempSTL && isfile(actualStlPath)
            delete(actualStlPath);
        end

        if wantOBJ && smoothed && colored
            info.obj = objPath;
            info.mtl = fullfile(fileparts(objPath), filenameStem + "_colored.mtl");
        end
    end

    % --- Voxel colored OBJ -----------------------------------------------
    if wantOBJ && ~smoothed
        effectiveLabels = labels;
        if isempty(effectiveLabels)
            effectiveLabels = zeros(numel(rho), 1);
        end
        exportDensityColoredOBJ(model, rho, effectiveLabels, labelNames, ...
            objPath, threshold, nModules, opts);
        info.obj = objPath;
        info.mtl = fullfile(fileparts(objPath), filenameStem + "_colored.mtl");
    end

    % --- MAT patch interchange -------------------------------------------
    if wantMAT
        exportCurveTopologyPatch(model, rho, matPath, threshold);
        info.mat = matPath;
    end
end

% =========================================================================
% Shared write helpers (used by 2D extrusion path)
% =========================================================================
function writeBinarySTL(filename, V, F)
    nTri    = size(F, 1);
    v1      = V(F(:,1),:);  v2 = V(F(:,2),:);  v3 = V(F(:,3),:);
    normals = cross(v2-v1, v3-v1, 2);
    normals = normals ./ max(sqrt(sum(normals.^2, 2)), eps);
    fid = fopen(filename, 'wb');
    if fid < 0
        error('exportTopology:cannotOpenSTL', 'Cannot open: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, sprintf('%-80s', 'Binary STL exported by exportTopology'), 'char');
    fwrite(fid, uint32(nTri), 'uint32');
    for i = 1:nTri
        fwrite(fid, single(normals(i,:)), 'single');
        fwrite(fid, single(v1(i,:)),      'single');
        fwrite(fid, single(v2(i,:)),      'single');
        fwrite(fid, single(v3(i,:)),      'single');
        fwrite(fid, uint16(0),            'uint16');
    end
end

function writeColoredOBJ(objPath, V, F, faceLabels, labelNames, colors, filenameStem)
    faceLabels = max(0, round(faceLabels(:)));
    names      = ["mixed"; string(labelNames(:))];
    maxLabel   = max([faceLabels; 0]);
    if maxLabel + 1 > numel(names)
        names = [names; "label" + string(numel(names):maxLabel)'];
    end
    if isempty(colors)
        colors = defaultOriginColors(numel(names));
    end
    mtlName = filenameStem + "_colored.mtl";
    mtlPath = fullfile(fileparts(objPath), mtlName);

    % Write MTL
    fid = fopen(mtlPath, 'w');
    if fid < 0, error('exportTopology:cannotOpenMTL', 'Cannot open: %s', mtlPath); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# MTL exported by exportTopology\n');
    for i = 1:numel(names)
        c = colors(min(i, size(colors,1)), :);
        mn = matName(names(i), i-1);
        fprintf(fid, 'newmtl %s\nKd %.6f %.6f %.6f\nKa %.6f %.6f %.6f\n', ...
            mn, c(1),c(2),c(3), 0.25*c(1),0.25*c(2),0.25*c(3));
        fprintf(fid, 'Ks 0.1 0.1 0.1\nNs 24\nd 1\n\n');
    end
    clear cleanup;

    % Write OBJ
    fid = fopen(objPath, 'w');
    if fid < 0, error('exportTopology:cannotOpenOBJ', 'Cannot open: %s', objPath); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# OBJ exported by exportTopology\nmtllib %s\n', mtlName);
    for i = 1:size(V,1)
        fprintf(fid, 'v %.9g %.9g %.9g\n', V(i,1), V(i,2), V(i,3));
    end
    fprintf(fid, '\n');
    for lab = unique(faceLabels)'
        mn  = matName(names(lab+1), lab);
        idx = find(faceLabels == lab);
        fprintf(fid, 'g %s\nusemtl %s\n', mn, mn);
        for k = 1:numel(idx)
            tri = F(idx(k),:);
            fprintf(fid, 'f %d %d %d\n', tri(1), tri(2), tri(3));
        end
        fprintf(fid, '\n');
    end
end

function mn = matName(name, label)
    raw = regexprep(char(name), '[^A-Za-z0-9_]+', '_');
    if isempty(raw), raw = sprintf('label%d', label); end
    mn = sprintf('mat_%02d_%s', label, raw);
end

function colors = defaultOriginColors(n)
    base = [
        0.55 0.55 0.55; 0.88 0.10 0.10; 0.10 0.25 0.90; 0.10 0.62 0.24
        0.95 0.55 0.05; 0.55 0.20 0.75; 0.95 0.85 0.05; 0.00 0.70 0.85
    ];
    if n <= size(base,1), colors = base(1:n,:);
    else, colors = [base; lines(n - size(base,1))];
    end
end

function tf = is2D(model)
    tf = size(model.mesh.nodes, 2) == 2;
end

function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name), value = s.(name); else, value = defaultValue; end
end

function info = exportTopology(model, rho, resultRoot, opts)
% exportTopology  Write topology geometry to external file formats.
%
%   exportTopology(model, rho, resultRoot)
%   exportTopology(model, rho, resultRoot, opts)
%   info = exportTopology(...)   % returns struct of written file paths
%
%   Supported formats
%     "stl"  Binary STL (smooth marching-tet or hard-threshold voxel)
%     "obj"  Colored OBJ + MTL (smooth or voxel, per-family colors)
%     "mat"  MATLAB patch data for later reimport
%
%   Format opts
%     formats          string | string array, subset of ["stl","obj","mat"]
%                      (default "stl")
%     smoothed         true | false — smooth marching-tet vs voxel (default true)
%     colored          true | false — per-family coloring for OBJ (default false)
%     threshold        iso-value (default 0.5)
%     nModules         restrict to first N modules (default: all)
%     filenameStem     base name for all written files (default 'topology')
%
%   Coloring opts (used when colored = true or labels is provided)
%     labels           [nElems x 1] integer label per element
%     labelNames       string array of label names
%     colors           [nLabels x 3] RGB
%
%   Smoothing opts (passed through to densitySmoothSurface when smoothed=true)
%     nodalSmoothingIters    (default 1)
%     surfaceSmoothingIters  (default 12)
%     taubinLambda           (default 0.45)
%     taubinMu               (default -0.47)
%
%   When both "stl" and "obj" are requested with smoothed=true, the smooth
%   geometry is computed once and the colored OBJ is written as a sidecar.

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

    % --- STL and smooth colored OBJ (share one geometry pass) ------------
    if wantSTL || (wantOBJ && smoothed)
        exportOpts              = opts;
        exportOpts.labels       = labels;
        exportOpts.labelNames   = labelNames;
        exportOpts.colors       = colors;

        % When OBJ is also requested, have the smooth exporter write it as
        % a sidecar so the geometry is only computed once.
        if wantOBJ && smoothed && colored
            exportOpts.coloredObjPath = objPath;
        end

        % If OBJ is the only smooth output (no STL), use a temp file.
        actualStlPath  = stlPath;
        deleteTempSTL  = false;
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
            info.mtl = fullfile(resultRoot, filenameStem + "_colored.mtl");
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
        info.mtl = fullfile(resultRoot, filenameStem + "_colored.mtl");
    end

    % --- MAT patch interchange -------------------------------------------
    if wantMAT
        exportCurveTopologyPatch(model, rho, matPath, threshold);
        info.mat = matPath;
    end
end

% -------------------------------------------------------------------------
function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name), value = s.(name); else, value = defaultValue; end
end

function h = plotDensitySmoothSurface(model, rho, opts)
% plotDensitySmoothSurface  Render a smooth density isosurface in the current axes.
%
%   h = plotDensitySmoothSurface(model, rho)
%   h = plotDensitySmoothSurface(model, rho, opts)
%
%   Extracts a marching-tetrahedra isosurface, applies Taubin smoothing,
%   and renders it as a patch object. Accepts the same smoothing options as
%   densitySmoothSurface plus the display options below.
%
%   Smoothing opts (passed through to densitySmoothSurface)
%     threshold             iso-value (default 0.5)
%     nModules              restrict to first N modules (default: all)
%     nodalSmoothingIters   (default 1)
%     surfaceSmoothingIters (default 12)
%     taubinLambda          (default 0.45)
%     taubinMu              (default -0.47)
%     labels                [nElems x 1] integer label per element
%     labelNames            string array of label names
%
%   Display opts
%     faceColor             [R G B] solid color when no labels (default [0.60 0.60 0.65])
%     colors                [nLabels x 3] RGB per label (default: built-in palette)
%     faceAlpha             opacity 0..1 (default 1.0)
%     edgeColor             edge color or 'none' (default 'none')
%     lighting              'gouraud' | 'flat' | 'none' (default 'gouraud')
%     specularStrength      (default 0.15)

    if nargin < 3 || isempty(opts), opts = struct(); end

    threshold = localOpt(opts, 'threshold', 0.5);
    nModules  = localOpt(opts, 'nModules',  0);
    faceAlpha = localOpt(opts, 'faceAlpha', 1.0);
    edgeColor = localOpt(opts, 'edgeColor', 'none');
    faceColor = localOpt(opts, 'faceColor', [0.60 0.60 0.65]);
    lightingMode      = localOpt(opts, 'lighting',         'gouraud');
    specularStrength  = localOpt(opts, 'specularStrength', 0.15);

    useLabels  = isfield(opts, 'labels') && ~isempty(opts.labels);
    colors     = localOpt(opts, 'colors', []);
    labelNames = string(localOpt(opts, 'labelNames', strings(0, 1)));

    [V, F, faceLabels] = densitySmoothSurface(model, rho, threshold, nModules, opts);

    if isempty(F)
        warning('plotDensitySmoothSurface:noSurface', ...
            'No isosurface found at threshold %.3g.', threshold);
        h = [];
        return;
    end

    if useLabels
        names    = ["mixed"; labelNames(:)];
        maxLabel = max(faceLabels);
        if maxLabel + 1 > numel(names)
            names = [names; "label" + string(numel(names):maxLabel)'];
        end
        if isempty(colors)
            colors = defaultOriginColors(numel(names));
        end
        faceRGB = colors(faceLabels + 1, :);
        h = patch('Vertices', V, 'Faces', F, ...
            'FaceVertexCData', faceRGB, 'FaceColor', 'flat', ...
            'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, ...
            'SpecularStrength', specularStrength);
    else
        h = patch('Vertices', V, 'Faces', F, ...
            'FaceColor', faceColor, ...
            'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, ...
            'SpecularStrength', specularStrength);
    end

    if ~strcmp(lightingMode, 'none')
        lighting(gca, lightingMode);
        if isempty(findobj(gcf, 'Type', 'light'))
            camlight('headlight');
        end
    end
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

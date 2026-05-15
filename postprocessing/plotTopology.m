function fig = plotTopology(model, rho, resultRoot, opts)
% plotTopology  Create and optionally save a topology figure (2D or 3D).
%
%   plotTopology(model, rho, resultRoot)
%   plotTopology(model, rho, resultRoot, opts)
%   fig = plotTopology(...)
%
%   Auto-detects 2D (quad mesh) vs 3D (hex mesh) from model.mesh.nodes.
%   Pass resultRoot = [] to render without saving.
%
%   Smoothing opts
%     smoothed         true | false (default true)
%     threshold        iso-value (default 0.5)
%     nModules         restrict to first N modules (default: all)
%     nodalSmoothingIters    (default 1)
%     surfaceSmoothingIters  (default 12)
%
%   Coloring opts
%     colored          true | false (default false)
%     labels           [nElems x 1] integer label per element
%     labelNames       string array of label names
%     colors           [nLabels x 3] RGB
%
%   Figure opts (3D)
%     viewAngle        [azimuth elevation] (default [45 35])
%
%   Save opts
%     figTitle         title string (default '')
%     filenameStem     base name for saved files (default 'topology')
%     saveFig          save .fig file (default false)
%     saveSVG          save .svg  — best for 2D (default false)
%     savePDF          save .pdf  — best for 2D (default false)
%     saveEPS          save .eps  — best for 2D (default false)
%     visible          keep figure window open (default false)

    if nargin < 4 || isempty(opts), opts = struct(); end

    smoothed     = localOpt(opts, 'smoothed',     true);
    colored      = localOpt(opts, 'colored',      false);
    threshold    = localOpt(opts, 'threshold',    0.5);
    nModules     = localOpt(opts, 'nModules',     0);
    viewAngle    = localOpt(opts, 'viewAngle',    [45 35]);
    figTitle     = string(localOpt(opts, 'figTitle',     ''));
    filenameStem = string(localOpt(opts, 'filenameStem', 'topology'));
    saveFig      = localOpt(opts, 'saveFig',  false);
    saveSVG      = localOpt(opts, 'saveSVG',  false);
    savePDF      = localOpt(opts, 'savePDF',  false);
    saveEPS      = localOpt(opts, 'saveEPS',  false);
    visible      = localOpt(opts, 'visible',  false);

    fig = figure('Visible', onOff(visible));
    ax  = axes('Parent', fig);
    hold(ax, 'on');

    if is2D(model)
        axis(ax, 'equal'); axis(ax, 'off');
        render2D(ax, model, rho, smoothed, colored, threshold, nModules, opts);
    else
        axis(ax, 'off'); daspect(ax, [1 1 1]);
        view(ax, viewAngle(1), viewAngle(2));
        render3D(ax, fig, model, rho, smoothed, colored, threshold, nModules, opts);
    end

    if strlength(figTitle) > 0
        title(ax, figTitle, 'Interpreter', 'none');
    end

    if ~isempty(resultRoot)
        if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end
        base = fullfile(resultRoot, filenameStem);
        saveas(fig, base + ".png");
        if saveFig, savefig(fig, base + ".fig");             end
        if saveSVG, print(fig, base + ".svg",  '-dsvg');    end
        if savePDF, print(fig, base + ".pdf",  '-dpdf');    end
        if saveEPS, print(fig, base + ".eps",  '-depsc');   end
    end

    if ~visible
        close(fig);
        if nargout == 0, clear fig; end
    end
end

% -------------------------------------------------------------------------
function render2D(ax, model, rho, smoothed, colored, threshold, nModules, opts)
    labels     = localOpt(opts, 'labels',     []);
    labelNames = string(localOpt(opts, 'labelNames', strings(0,1)));
    colors     = localOpt(opts, 'colors',     []);

    if smoothed
        smoothOpts           = opts;
        smoothOpts.threshold = threshold;
        smoothOpts.nModules  = nModules;
        if ~colored, smoothOpts = rmfield_safe(smoothOpts, 'labels'); end

        [polygons, polyLabels] = densitySmoothSurface2D( ...
            model, rho, threshold, nModules, smoothOpts);

        if colored && ~isempty(polygons)
            names  = ["mixed"; labelNames(:)];
            nNames = max(max(polyLabels) + 2, numel(names));
            if isempty(colors), colors = defaultOriginColors(nNames); end
            for i = 1:numel(polygons)
                c = colors(polyLabels(i) + 1, :);
                fill(ax, polygons{i}(:,1), polygons{i}(:,2), c, 'EdgeColor', 'none');
            end
        else
            for i = 1:numel(polygons)
                fill(ax, polygons{i}(:,1), polygons{i}(:,2), ...
                    [0.60 0.60 0.65], 'EdgeColor', 'none');
            end
        end
    else
        [V, F, faceLabels] = densityVoxelSurface2D( ...
            model, rho, threshold, nModules, labels);
        if isempty(F), return; end

        if colored && ~isempty(labels)
            names  = ["mixed"; labelNames(:)];
            nNames = max(max(faceLabels) + 2, numel(names));
            if isempty(colors), colors = defaultOriginColors(nNames); end
            faceRGB = colors(faceLabels + 1, :);
            patch(ax, 'Vertices', V, 'Faces', F, ...
                'FaceVertexCData', faceRGB, 'FaceColor', 'flat', 'EdgeColor', 'none');
        else
            patch(ax, 'Vertices', V, 'Faces', F, ...
                'FaceColor', [0.60 0.60 0.65], 'EdgeColor', 'none');
        end
    end
end

% -------------------------------------------------------------------------
function render3D(ax, fig, model, rho, smoothed, colored, threshold, nModules, opts)
    if smoothed
        renderOpts           = opts;
        renderOpts.threshold = threshold;
        renderOpts.nModules  = nModules;
        if ~colored, renderOpts = rmfield_safe(renderOpts, 'labels'); end
        plotDensitySmoothSurface(model, rho, renderOpts);
    else
        if colored
            warning('plotTopology:coloredVoxelNotSupported', ...
                'colored=true requires smoothed=true for 3D; rendering gray voxel surface.');
        end
        [V, F] = densityVoxelSurface(model, rho, threshold, nModules);
        if ~isempty(F)
            patch(ax, 'Vertices', V, 'Faces', F, ...
                'FaceColor', [0.60 0.60 0.65], 'EdgeColor', 'none');
            lighting(ax, 'gouraud');
            if isempty(findobj(fig, 'Type', 'light'))
                camlight(ax, 'headlight');
            end
        end
    end
end

% -------------------------------------------------------------------------
function tf = is2D(model)
    tf = size(model.mesh.nodes, 2) == 2;
end

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

function s = rmfield_safe(s, name)
    if isfield(s, name), s = rmfield(s, name); end
end

function str = onOff(flag)
    if flag, str = 'on'; else, str = 'off'; end
end

function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name), value = s.(name); else, value = defaultValue; end
end

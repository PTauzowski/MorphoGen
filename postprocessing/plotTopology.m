function fig = plotTopology(model, rho, resultRoot, opts)
% plotTopology  Create and optionally save a MATLAB topology figure.
%
%   plotTopology(model, rho, resultRoot)
%   plotTopology(model, rho, resultRoot, opts)
%   fig = plotTopology(...)
%
%   Renders the density field as a surface patch and saves PNG (and
%   optionally FIG) to resultRoot. Pass resultRoot = [] to skip saving.
%
%   Smoothing opts (passed through to densitySmoothSurface when smoothed=true)
%     smoothed         true | false (default true)
%     threshold        iso-value (default 0.5)
%     nModules         restrict to first N modules (default: all)
%     nodalSmoothingIters    (default 1)
%     surfaceSmoothingIters  (default 12)
%
%   Coloring opts (only effective when smoothed = true)
%     colored          true | false (default false)
%     labels           [nElems x 1] integer label per element
%     labelNames       string array of label names
%     colors           [nLabels x 3] RGB
%
%   Figure opts
%     viewAngle        [azimuth elevation] (default [45 35])
%     figTitle         title string (default '')
%     filenameStem     base name for saved files (default 'topology')
%     saveFig          also save .fig file (default false)
%     visible          keep figure window open (default false)

    if nargin < 4 || isempty(opts), opts = struct(); end

    smoothed     = localOpt(opts, 'smoothed',     true);
    colored      = localOpt(opts, 'colored',      false);
    threshold    = localOpt(opts, 'threshold',    0.5);
    nModules     = localOpt(opts, 'nModules',     0);
    viewAngle    = localOpt(opts, 'viewAngle',    [45 35]);
    figTitle     = string(localOpt(opts, 'figTitle',     ''));
    filenameStem = string(localOpt(opts, 'filenameStem', 'topology'));
    saveFig      = localOpt(opts, 'saveFig',      false);
    visible      = localOpt(opts, 'visible',      false);

    fig = figure('Visible', onOff(visible));
    ax  = axes('Parent', fig);
    hold(ax, 'on'); axis(ax, 'off'); daspect(ax, [1 1 1]);
    view(ax, viewAngle(1), viewAngle(2));

    if smoothed
        renderOpts           = opts;
        renderOpts.threshold = threshold;
        renderOpts.nModules  = nModules;
        if ~colored
            renderOpts = rmfield_safe(renderOpts, 'labels');
        end
        plotDensitySmoothSurface(model, rho, renderOpts);
    else
        if colored
            warning('plotTopology:coloredVoxelNotSupported', ...
                'colored=true requires smoothed=true; rendering gray voxel surface.');
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

    if strlength(figTitle) > 0
        title(ax, figTitle, 'Interpreter', 'none');
    end

    if ~isempty(resultRoot)
        if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end
        saveas(fig, fullfile(resultRoot, filenameStem + ".png"));
        if saveFig
            savefig(fig, fullfile(resultRoot, filenameStem + ".fig"));
        end
    end

    if ~visible
        close(fig);
        if nargout == 0, clear fig; end
    end
end

% -------------------------------------------------------------------------
function s = rmfield_safe(s, name)
    if isfield(s, name), s = rmfield(s, name); end
end

function str = onOff(flag)
    if flag, str = 'on'; else, str = 'off'; end
end

function value = localOpt(s, name, defaultValue)
    if isstruct(s) && isfield(s, name), value = s.(name); else, value = defaultValue; end
end

function fig = plotLinkedElement(model, varargin)
% PLOTLINKEDELEMENT  Prove the segmentToArm linkage by colouring all
%   full-arm copies of one reference half-segment element.
%
%   fig = plotLinkedElement(model)
%       Picks a random element e_ref from [1..H].
%
%   fig = plotLinkedElement(model, e_ref)
%       Uses the specified reference element index.
%
%   fig = plotLinkedElement(model, e_ref, resultRoot)
%       Also saves the figure to resultRoot as
%       'linked_element_<e_ref>.png' and '.fig'.
%
%   fig = plotLinkedElement(model, e_ref, resultRoot, configLabel, betas)
%       Adds the configuration label and beta sequence to the figure title.
%
%   INPUTS
%     model      ManipulatorModel3D  Full-arm solid model.
%     e_ref      (optional) Reference half-segment element index in [1..H].
%                Default: randomly chosen.
%     resultRoot (optional) Directory path for saving output files.
%
%   OUTPUT
%     fig        Figure handle (visible on screen).
%
%   WHAT IT SHOWS
%     The full arm mesh is drawn as a light-grey ghost.  The linked
%     elements — those that share the same design variable rho(e_ref)
%     via the segmentToArm expansion — are highlighted in two colours:
%       RED   : normal copies  (odd half-segment blocks, 2a-type)
%       BLUE  : flipped copies (even half-segment blocks, 2b-type,
%               the mirror partner at position H+1-e_ref)
%     Each segment boundary ring is shown as a thin dark line to help
%     count the copies visually.

    % ---- parse inputs -------------------------------------------------------
    e_ref      = [];
    resultRoot = '';
    configLabel = "";
    betas = [];
    if nargin >= 2 && ~isempty(varargin{1})
        e_ref = varargin{1};
    end
    if nargin >= 3
        resultRoot = varargin{2};
    end
    if nargin >= 4 && ~isempty(varargin{3})
        configLabel = string(varargin{3});
    end
    if nargin >= 5 && ~isempty(varargin{4})
        betas = varargin{4};
    end

    H      = model.halfSegmentNelems;
    nElems = size(model.mesh.elems, 1);
    assert(mod(nElems, H) == 0, ...
        'Full-arm element count %d is not divisible by H=%d.', nElems, H);
    nCopies = nElems / H;
    assert(mod(nCopies, 2) == 0, ...
        'nCopies = %d must be even (half-segment pairs).', nCopies);
    nArms = nCopies / 2;

    if isempty(e_ref)
        rng('shuffle');
        e_ref = randi(H);
    end
    assert(e_ref >= 1 && e_ref <= H, ...
        'e_ref=%d is out of range [1..%d].', e_ref, H);

    % ---- build element index sets via the sensitivity map -------------------
    map = buildLinkedSensitivityMap(H, nElems);

    % Normal copies: map.normalMap(:, e_ref)   — one per arm, column vector
    % Flipped copies: map.flippedMap(:, e_ref) — one per arm, column vector
    normalIds  = map.normalMap(:, e_ref);    % [nArms x 1]
    flippedIds = map.flippedMap(:, e_ref);   % [nArms x 1]

    % Boolean masks over the full element list
    maskNormal  = false(nElems, 1);
    maskFlipped = false(nElems, 1);
    maskNormal(normalIds)  = true;
    maskFlipped(flippedIds) = true;
    normalCentroids = elementCentroids(model.mesh.nodes, model.mesh.elems, normalIds);
    flippedCentroids = elementCentroids(model.mesh.nodes, model.mesh.elems, flippedIds);
    configText = "";
    if strlength(configLabel) > 0
        configText = configLabel;
    end
    if ~isempty(betas)
        betaText = "betas=" + string(mat2str(betas));
        if strlength(configText) > 0
            configText = configText + " | " + betaText;
        else
            configText = betaText;
        end
    end

    % ---- save fe rendering state so we can restore it afterwards -----------
    savedAlpha = model.fe.face_alpha;
    savedFaceColor = model.fe.face_color;
    savedEdgeColor = model.fe.edge_color;

    % ---- build figure -------------------------------------------------------
    fig = figure('Name', sprintf('Linked element e_ref=%d (H=%d, nArms=%d)', ...
        e_ref, H, nArms), ...
        'Color', 'white', 'Units', 'normalized', 'Position', [0.05 0.05 0.6 0.85]);

    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ---- Tile 1 : 3-D perspective -------------------------------------------
    nexttile;
    hold on; axis off; daspect([1 1 1]);
    view(45, 30);
    light('Position', [-1 -1 3], 'Style', 'infinite');
    light('Position', [1  2  1], 'Style', 'local');

    % Ghost: full arm, very transparent grey
    model.fe.face_alpha = 0.08;
    model.fe.face_color = [0.6 0.6 0.6];
    model.fe.edge_color = 'none';
    model.fe.plot(model.mesh.nodes);

    % Normal copies — red
    model.fe.face_alpha = 1.0;
    model.fe.face_color = [0.85 0.15 0.10];
    model.fe.edge_color = [0.5  0.05 0.05];
    model.fe.plotSolidSelected(model.mesh.nodes, maskNormal, [0.85 0.15 0.10]);

    % Flipped copies — blue
    model.fe.face_color = [0.15 0.35 0.85];
    model.fe.edge_color = [0.05 0.10 0.50];
    model.fe.plotSolidSelected(model.mesh.nodes, maskFlipped, [0.15 0.35 0.85]);

    model.fe.face_alpha = 1.0;
    model.fe.edge_color = [0.3 0.3 0.3];
    addElementLabels(normalCentroids, "N", [0.55 0.02 0.02]);
    addElementLabels(flippedCentroids, "F", [0.02 0.08 0.50]);

    if strlength(configText) > 0
        title(sprintf('3-D view  |  %s  |  e_{ref}=%d, H=%d, %d copies', ...
            configText, e_ref, H, nCopies), ...
            'Interpreter', 'none');
    else
        title(sprintf('3-D view  (e_{ref}=%d, H=%d, %d copies)', e_ref, H, nCopies), ...
        'Interpreter', 'tex');
    end
    xlabel('x'); ylabel('y'); zlabel('z');

    % ---- Tile 2 : side view (YZ plane) --------------------------------------
    nexttile;
    hold on; axis off; daspect([1 1 1]);
    view(90, 0);
    light('Position', [0 -1 1], 'Style', 'infinite');

    model.fe.face_alpha = 0.08;
    model.fe.face_color = [0.6 0.6 0.6];
    model.fe.edge_color = 'none';
    model.fe.plot(model.mesh.nodes);

    model.fe.face_alpha = 1.0;
    model.fe.face_color = [0.85 0.15 0.10];
    model.fe.edge_color = [0.5  0.05 0.05];
    model.fe.plotSolidSelected(model.mesh.nodes, maskNormal, [0.85 0.15 0.10]);

    model.fe.face_color = [0.15 0.35 0.85];
    model.fe.edge_color = [0.05 0.10 0.50];
    model.fe.plotSolidSelected(model.mesh.nodes, maskFlipped, [0.15 0.35 0.85]);

    model.fe.face_alpha = 1.0;
    model.fe.edge_color = [0.3 0.3 0.3];
    addElementLabels(normalCentroids, "N", [0.55 0.02 0.02]);
    addElementLabels(flippedCentroids, "F", [0.02 0.08 0.50]);

    if strlength(configText) > 0
        title(sprintf('Side view (YZ)  |  %s', configText), 'Interpreter', 'none');
    else
        title('Side view (YZ)', 'Interpreter', 'tex');
    end
    xlabel('x'); ylabel('y'); zlabel('z');

    % ---- super-title --------------------------------------------------------
    mirrorId = H + 1 - e_ref;
    if strlength(configText) > 0
        sgtitle(sprintf( ...
            ['Linked element proof  |  %s\n' ...
             'rho(%d) drives %d red + %d blue = %d elements  |  Red=N normal copies, Blue=F flipped copies'], ...
            configText, e_ref, nArms, nArms, nCopies), ...
            'Interpreter', 'none', 'FontSize', 11);
    else
        sgtitle(sprintf( ...
            ['Linked element proof  |  \\rho(%d)  drives  %d red + %d blue = %d elements\n' ...
             'Red = e_{ref}=%d in normal (2a) copies    Blue = e_{ref}=%d in flipped (2b) copies'], ...
            e_ref, nArms, nArms, nCopies, e_ref, mirrorId), ...
            'Interpreter', 'tex', 'FontSize', 11);
    end

    % ---- print stats to console ---------------------------------------------
    fprintf('\nplotLinkedElement: e_ref=%d  (H=%d, nCopies=%d, nArms=%d)\n', ...
        e_ref, H, nCopies, nArms);
    fprintf('  Normal  copies (2a, red)  : %d elements at global ids: %s\n', ...
        numel(normalIds),  mat2str(normalIds(:)'));
    fprintf('  Flipped copies (2b, blue) : %d elements at global ids: %s\n', ...
        numel(flippedIds), mat2str(flippedIds(:)'));
    fprintf('  Mirror element in rho-space: e_ref=%d  <->  H+1-e_ref=%d\n', ...
        e_ref, mirrorId);
    if strlength(configText) > 0
        fprintf('  Configuration: %s\n', configText);
    end

    % ---- restore fe rendering state -----------------------------------------
    model.fe.face_alpha = savedAlpha;
    model.fe.face_color = savedFaceColor;
    model.fe.edge_color = savedEdgeColor;

    % ---- optional save ------------------------------------------------------
    if ~isempty(resultRoot)
        if ~exist(resultRoot, 'dir')
            mkdir(resultRoot);
        end
        stem = fullfile(resultRoot, sprintf('linked_element_%d', e_ref));
        exportgraphics(fig, [stem '.png'], 'Resolution', 200);
        savefig(fig, [stem '.fig']);
        fprintf('  Saved to %s\n', stem);
    end
end

function centroids = elementCentroids(nodes, elems, elemIds)
    centroids = zeros(numel(elemIds), size(nodes, 2));
    for i = 1:numel(elemIds)
        centroids(i, :) = mean(nodes(elems(elemIds(i), :), :), 1);
    end
end

function addElementLabels(centroids, prefix, color)
    for i = 1:size(centroids, 1)
        text(centroids(i,1), centroids(i,2), centroids(i,3), ...
            sprintf('%s%d', prefix, i), ...
            'Color', color, ...
            'FontSize', 8, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'BackgroundColor', 'white', ...
            'Margin', 1, ...
            'Clipping', 'on');
    end
end

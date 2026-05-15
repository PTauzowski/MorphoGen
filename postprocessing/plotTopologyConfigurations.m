function plotTopologyConfigurations(models, elemMask, resultRoot, filenameStem, titleText, configs)
% plotTopologyConfigurations  Plot one full-arm element mask on each config.
%
% Linked tests use one element-density vector for all configurations. The
% ordinary postprocess figures show only models{1}; this diagnostic renders
% the same mask on every configuration so local-frame rotations are visible.

    if nargin < 6
        configs = [];
    end

    filenameStem = string(filenameStem);
    saveFigFiles = false;
    nConfigs = numel(models);
    elemMask = logical(elemMask(:));
    assert(nConfigs > 0, 'plotTopologyConfigurations: models must be non-empty.');
    assert(numel(elemMask) == models{1}.analysis.getTotalElemsNumber(), ...
        'plotTopologyConfigurations: mask length %d does not match element count %d.', ...
        numel(elemMask), models{1}.analysis.getTotalElemsNumber());

    fig = figure('Visible', 'off', 'Name', filenameStem);
    nCols = min(3, nConfigs);
    nRows = ceil(nConfigs / nCols);
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:nConfigs
        nexttile;
        hold on; axis off; daspect([1 1 1]); view(45, 35);
        models{k}.fe.plotSolidSelected(models{k}.mesh.nodes, elemMask, [0.45 0.60 0.80]);

        if ~isempty(configs)
            label = string(configs{k}.label);
        else
            label = "config " + k;
        end
        title(label, 'Interpreter', 'none');
    end

    sgtitle(titleText, 'Interpreter', 'none');
    saveas(fig, fullfile(resultRoot, filenameStem + ".png"));
    if saveFigFiles
        savefig(fig, fullfile(resultRoot, filenameStem + ".fig"));
    end
    close(fig);
end

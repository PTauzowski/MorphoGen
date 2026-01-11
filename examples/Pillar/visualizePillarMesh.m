function visualizePillarMesh(mesh, opts)
% VISUALIZEPILLARMESH - Create visualization plots for pillar mesh
%
% SYNTAX:
%   visualizePillarMesh(mesh)
%   visualizePillarMesh(mesh, opts)
%
% INPUTS:
%   mesh - Mesh object from generatePillarQuarterMesh
%   opts - Optional visualization options:
%     showMaterials - Color by material ID (default: true)
%     viewAngle - 3D view angle [az, el] (default: [45, 30])
%     sliceZ - z-level for horizontal slice plot (default: 0)
%
% DESCRIPTION:
%   Creates a comprehensive visualization of the pillar mesh including:
%   - 3D view of all nodes
%   - Horizontal slice at specified z
%   - Vertical slice (x-z plane)
%   - Material distribution
%   - Node statistics

    if nargin < 2
        opts = struct();
    end

    % Set defaults
    if ~isfield(opts, 'showMaterials')
        opts.showMaterials = true;
    end
    if ~isfield(opts, 'viewAngle')
        opts.viewAngle = [45, 30];
    end
    if ~isfield(opts, 'sliceZ')
        opts.sliceZ = 0;
    end

    % Extract data
    nodes = mesh.nodes;
    elems = mesh.elems;

    % Create figure
    figure('Name', 'Pillar Mesh Visualization', 'Position', [50, 50, 1400, 900]);

    % Plot 1: 3D view colored by z
    subplot(2, 3, 1);
    scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 2, nodes(:,3), 'filled');
    axis equal; grid on; box on;
    xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
    title('3D View (colored by z)');
    view(opts.viewAngle);
    colorbar;

    % Plot 2: Top view (x-y plane, all nodes)
    subplot(2, 3, 2);
    scatter(nodes(:,1), nodes(:,2), 3, nodes(:,3), 'filled');
    axis equal; grid on; box on;
    xlabel('x [nm]'); ylabel('y [nm]');
    title('Top View (colored by z)');
    colorbar;

    % Plot 3: Side view (x-z plane, near y=0)
    subplot(2, 3, 3);
    tol = max(nodes(:,2)) * 0.05; % 5% tolerance
    idx_y0 = abs(nodes(:,2)) < tol;
    scatter(nodes(idx_y0,1), nodes(idx_y0,3), 4, 'b', 'filled');
    grid on; box on;
    xlabel('x [nm]'); ylabel('z [nm]');
    title(sprintf('Side View (|y| < %.1f nm)', tol));

    % Plot 4: Horizontal slice at specified z
    subplot(2, 3, 4);
    tol_z = max(abs(nodes(:,3))) * 0.01; % 1% tolerance
    idx_slice = abs(nodes(:,3) - opts.sliceZ) < tol_z;
    if sum(idx_slice) > 0
        scatter(nodes(idx_slice,1), nodes(idx_slice,2), 20, 'r', 'filled');
        axis equal; grid on; box on;
        xlabel('x [nm]'); ylabel('y [nm]');
        title(sprintf('Horizontal Slice (z ≈ %.1f nm)', opts.sliceZ));
    else
        text(0.5, 0.5, 'No nodes at this z-level', ...
            'HorizontalAlignment', 'center');
        title(sprintf('Horizontal Slice (z ≈ %.1f nm) - Empty', opts.sliceZ));
    end

    % Plot 5: Material distribution
    subplot(2, 3, 5);
    if isfield(mesh, 'matID') && ~isempty(mesh.matID)
        % Compute element centroids
        nElems = size(elems, 1);
        centroids = zeros(nElems, 3);
        for i = 1:nElems
            centroids(i, :) = mean(nodes(elems(i, :), :), 1);
        end

        % Plot colored by material
        scatter3(centroids(:,1), centroids(:,2), centroids(:,3), ...
            10, mesh.matID, 'filled');
        axis equal; grid on; box on;
        xlabel('x [nm]'); ylabel('y [nm]'); zlabel('z [nm]');
        title('Element Centroids (colored by material ID)');
        view(opts.viewAngle);
        colorbar;

        % Add material ID legend info
        uniqueMats = unique(mesh.matID);
        legendText = sprintf('Material IDs: %s', mat2str(uniqueMats'));
        text(0.02, 0.98, legendText, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'BackgroundColor', 'w');
    else
        text(0.5, 0.5, 'No material IDs available', ...
            'HorizontalAlignment', 'center');
        title('Material Distribution - Not Available');
    end

    % Plot 6: Statistics and layer information
    subplot(2, 3, 6);
    axis off;

    % Compute statistics
    stats = {};
    stats{end+1} = sprintf('MESH STATISTICS');
    stats{end+1} = sprintf('====================');
    stats{end+1} = sprintf('Nodes: %d', size(nodes, 1));
    stats{end+1} = sprintf('Elements: %d', size(elems, 1));
    stats{end+1} = sprintf('Nodes per element: %d', size(elems, 2));
    stats{end+1} = '';

    if isfield(mesh, 'sf')
        stats{end+1} = sprintf('Shape function: %s', class(mesh.sf));
    end
    stats{end+1} = '';

    stats{end+1} = sprintf('COORDINATE RANGES');
    stats{end+1} = sprintf('====================');
    stats{end+1} = sprintf('x: [%.2f, %.2f] nm', min(nodes(:,1)), max(nodes(:,1)));
    stats{end+1} = sprintf('y: [%.2f, %.2f] nm', min(nodes(:,2)), max(nodes(:,2)));
    stats{end+1} = sprintf('z: [%.2f, %.2f] nm', min(nodes(:,3)), max(nodes(:,3)));
    stats{end+1} = '';

    % Z-distribution
    uniqueZ = unique(nodes(:,3));
    stats{end+1} = sprintf('Z-LEVELS');
    stats{end+1} = sprintf('====================');
    stats{end+1} = sprintf('Unique z-values: %d', length(uniqueZ));
    stats{end+1} = sprintf('Min spacing: %.3f nm', min(diff(sort(uniqueZ))));
    stats{end+1} = sprintf('Max spacing: %.3f nm', max(diff(sort(uniqueZ))));
    stats{end+1} = '';

    % Material info
    if isfield(mesh, 'matID') && ~isempty(mesh.matID)
        stats{end+1} = sprintf('MATERIALS');
        stats{end+1} = sprintf('====================');
        uniqueMats = unique(mesh.matID);
        for i = 1:length(uniqueMats)
            matID = uniqueMats(i);
            nElems = sum(mesh.matID == matID);
            stats{end+1} = sprintf('ID %d: %d elements', matID, nElems);
        end
    end

    % Display text
    text(0.05, 0.95, strjoin(stats, '\n'), ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'FixedWidth', ...
        'FontSize', 8);

    % Overall title
    sgtitle('Pillar Quarter Mesh Visualization', 'FontSize', 14, 'FontWeight', 'bold');
end

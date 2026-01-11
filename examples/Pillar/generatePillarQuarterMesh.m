function mesh = generatePillarQuarterMesh(opts)
% GENERATEPILLARQUARTERMESH - Generate 3D hexahedral mesh for quarter axisymmetric pillar-on-ground
%
% SYNTAX:
%   mesh = generatePillarQuarterMesh(opts)
%
% DESCRIPTION:
%   Generates a quarter-tile mesh of a pillar structure with:
%   - Cylindrical-like mesh near pillar (core + bank)
%   - Transition region warping circle→square toward outer boundary
%   - Proper layer stack with horizontal interfaces (NO z-warping)
%   - Etched regions represented as void (no elements) or optional air material
%
% PARAMETERS (opts structure with defaults):
%   Pillar geometry:
%     rTop          - Pillar top radius [nm] (default: 100)
%     alpha_deg     - Cone inclination from vertical [deg] (default: 11)
%     zBank         - Bank/cap height above z=0 [nm] (default: 90)
%
%   Domain size:
%     Rout          - Outer domain radius [nm] (default: 1000)
%     zSubExtra     - Extra substrate below -350 nm (default: 0)
%
%   Trench/bank profile:
%     Rtrench       - Outer radius of trench zone [nm] (default: 300)
%     trenchWidth   - Radial width of trench [nm] (default: 150)
%     trenchDepth   - Maximum trench depth below z=0 [nm] (default: 200)
%     Hcurve        - Vertical smoothing height for bank [nm] (default: 50)
%     trenchSharpness - Sharpness of trench profile (default: 2.0)
%
%   Interface sublayers:
%     tInterface    - Interface layer thickness [nm] (default: 0.1)
%
%   Mesh resolution:
%     target_dz     - Target element size in z [nm] (default: 20)
%     nr_core       - Number of radial divisions in core (default: 8)
%     nr_trench     - Number of radial divisions in trench (default: 8)
%     nr_outer      - Number of radial divisions in outer zone (default: 6)
%     ntheta        - Number of angular divisions in quarter (default: 8)
%
%   Element type:
%     sfName        - Shape function name: 'H27' or 'H8' (default: 'H27')
%
%   Etched region mode:
%     etch_mode     - 'void' (no elements) or 'air_material' (default: 'void')
%
% OUTPUT:
%   mesh - Structure with fields:
%     mesh.nodes     - (N x 3) node coordinates [x, y, z]
%     mesh.elems     - (E x nen) element connectivity
%     mesh.matID     - (E x 1) material IDs for each element
%     mesh.sf        - ShapeFunction object used
%     mesh.meshObj   - Original Mesh object (for use with Mesh methods)
%
% PHYSICAL INTERPRETATION:
%   1) Flat horizontal layer stack deposited on substrate
%   2) Pillar and recess formed by etching (modifies solid/void boundary only)
%   3) All layer interfaces remain horizontal planes at fixed z-levels
%   4) z=0 is top of porous layer
%   5) Pillar is cone with 11° inclination, expanding downward from top
%
% EXAMPLE:
%   % Generate with defaults
%   mesh = generatePillarQuarterMesh();
%
%   % Custom resolution
%   opts.target_dz = 15;
%   opts.nr_core = 12;
%   mesh = generatePillarQuarterMesh(opts);
%
% AUTHOR: Generated for MorphoGenVibrations framework

    % Set default parameters
    if nargin < 1, opts = struct(); end

    % Pillar geometry
    opts = setDefault(opts, 'rTop', 100);
    opts = setDefault(opts, 'alpha_deg', 11);
    opts = setDefault(opts, 'zBank', 90);

    % Domain
    opts = setDefault(opts, 'Rout', 1000);
    opts = setDefault(opts, 'zSubExtra', 0);

    % Trench/bank
    opts = setDefault(opts, 'Rtrench', 300);
    opts = setDefault(opts, 'trenchWidth', 150);
    opts = setDefault(opts, 'trenchDepth', 200);
    opts = setDefault(opts, 'Hcurve', 50);
    opts = setDefault(opts, 'trenchSharpness', 2.0);

    % Interfaces
    opts = setDefault(opts, 'tInterface', 0.1);

    % Mesh resolution
    opts = setDefault(opts, 'target_dz', 20);
    opts = setDefault(opts, 'nr_core', 8);
    opts = setDefault(opts, 'nr_trench', 8);
    opts = setDefault(opts, 'nr_outer', 6);
    opts = setDefault(opts, 'ntheta', 8);

    % Element type
    opts = setDefault(opts, 'sfName', 'H27');

    % Etched region
    opts = setDefault(opts, 'etch_mode', 'void');

    % Build the mesh
    mesh = buildPillarMesh(opts);
end

function opts = setDefault(opts, field, value)
    if ~isfield(opts, field)
        opts.(field) = value;
    end
end

function mesh = buildPillarMesh(opts)
    % Main mesh builder

    fprintf('=== Pillar Mesh Generator ===\n');
    fprintf('Configuration:\n');
    fprintf('  Pillar top radius: %.1f nm\n', opts.rTop);
    fprintf('  Cone angle: %.1f deg\n', opts.alpha_deg);
    fprintf('  Bank height: %.1f nm\n', opts.zBank);
    fprintf('  Domain size: %.1f nm\n', opts.Rout);
    fprintf('  Shape function: %s\n', opts.sfName);
    fprintf('\n');

    % Build vertical layer stack (authoritative z-levels)
    [layers, zLevels, zTop] = buildLayerStack(opts);

    fprintf('Layer stack built:\n');
    fprintf('  Physical layers: %d\n', sum([layers.isInterface] == false));
    fprintf('  Interface layers: %d\n', sum([layers.isInterface] == true));
    fprintf('  Total z-levels: %d\n', length(zLevels));
    fprintf('  Pillar top z: %.2f nm\n', zTop);
    if any(abs(zLevels) < 1e-10)
        fprintf('  Verify z=0 in stack: true\n');
    else
        fprintf('  Verify z=0 in stack: false\n');
    end
    if any(abs(zLevels - opts.zBank) < 1e-10)
        fprintf('  Verify zBank in stack: true\n');
    else
        fprintf('  Verify zBank in stack: false\n');
    end
    fprintf('\n');

    % Create shape function
    sf = createShapeFunction(opts.sfName);

    % Define radial zones
    alpha_rad = opts.alpha_deg * pi / 180;
    r_pillar_at_z0 = opts.rTop + zTop * tan(alpha_rad);
    Rp_zone = r_pillar_at_z0 * 1.2; % Core zone slightly larger than pillar base

    fprintf('Radial zones:\n');
    fprintf('  Core: 0 -> %.1f nm\n', Rp_zone);
    fprintf('  Trench: %.1f -> %.1f nm\n', Rp_zone, opts.Rtrench);
    fprintf('  Outer: %.1f -> %.1f nm\n', opts.Rtrench, opts.Rout);
    fprintf('\n');

    % Generate mesh blocks
    fprintf('Generating mesh blocks...\n');
    blocks = {};

    % Core cylindrical zone
    fprintf('  Core zone...\n');
    blockCore = generateCylindricalBlock(0, Rp_zone, zLevels, ...
        opts.nr_core, opts.ntheta, layers, opts, sf, zTop);
    if ~isempty(blockCore.nodes)
        blocks{end+1} = blockCore;
    end

    % Trench annulus
    fprintf('  Trench zone...\n');
    blockTrench = generateCylindricalBlock(Rp_zone, opts.Rtrench, zLevels, ...
        opts.nr_trench, opts.ntheta, layers, opts, sf, zTop);
    if ~isempty(blockTrench.nodes)
        blocks{end+1} = blockTrench;
    end

    % Outer transition zone (circle to square)
    fprintf('  Outer transition zone...\n');
    blockOuter = generateTransitionBlock(opts.Rtrench, opts.Rout, zLevels, ...
        opts.nr_outer, opts.ntheta, layers, opts, sf, zTop);
    if ~isempty(blockOuter.nodes)
        blocks{end+1} = blockOuter;
    end

    % Merge all blocks
    fprintf('Merging blocks...\n');
    meshObj = Mesh();
    allMatID = [];

    for i = 1:length(blocks)
        meshObj.merge(blocks{i}.nodes, blocks{i}.elems);
        allMatID = [allMatID; blocks{i}.matID];
    end

    fprintf('Mesh generation complete.\n');
    fprintf('  Total nodes: %d\n', size(meshObj.nodes, 1));
    fprintf('  Total elements: %d\n', size(meshObj.elems, 1));
    fprintf('\n');

    % Create output structure with additional properties
    % We can't add properties to Mesh class, so create a wrapper structure
    mesh = struct();
    mesh.nodes = meshObj.nodes;
    mesh.elems = meshObj.elems;
    mesh.matID = allMatID;
    mesh.sf = sf;
    mesh.meshObj = meshObj;  % Keep reference to original Mesh object

    % Verify mesh
    verifyMesh(mesh, opts, zTop, zLevels, layers);
end

function [layers, zLevels, zTop] = buildLayerStack(opts)
    % Build authoritative vertical layer stack with interfaces
    % All layers are horizontal; z=0 and zBank are exact slab boundaries

    % Physical layer definitions (bottom to top)
    % Below z=0:
    physLayers = struct([]);
    layerIdx = 0;
    z_current = -350 - opts.zSubExtra;

    if opts.zSubExtra > 0
        layerIdx = layerIdx + 1;
        physLayers(layerIdx).name = 'Substrate_extra';
        physLayers(layerIdx).z0 = z_current;
        physLayers(layerIdx).thickness = opts.zSubExtra;
        physLayers(layerIdx).matID = 1;
        physLayers(layerIdx).nz = max(1, round(opts.zSubExtra / opts.target_dz));
        z_current = z_current + opts.zSubExtra;
    end

    % 50 nm GaN (non-porous)
    layerIdx = layerIdx + 1;
    physLayers(layerIdx).name = 'GaN_substrate';
    physLayers(layerIdx).z0 = z_current;
    physLayers(layerIdx).thickness = 50;
    physLayers(layerIdx).matID = 2;
    physLayers(layerIdx).nz = max(1, round(50 / opts.target_dz));
    z_current = z_current + 50;

    % 300 nm porous GaN (ends exactly at z=0)
    layerIdx = layerIdx + 1;
    physLayers(layerIdx).name = 'GaN_porous';
    physLayers(layerIdx).z0 = z_current;
    physLayers(layerIdx).thickness = 300;
    physLayers(layerIdx).matID = 3;
    physLayers(layerIdx).nz = max(2, round(300 / opts.target_dz));
    z_current = z_current + 300;

    % Verify we're at z=0
    assert(abs(z_current) < 1e-10, 'Stack must reach exactly z=0');

    % Above z=0 (pillar stack):
    % 10 nm GaN
    layerIdx = layerIdx + 1;
    physLayers(layerIdx).name = 'GaN_10nm';
    physLayers(layerIdx).z0 = z_current;
    physLayers(layerIdx).thickness = 10;
    physLayers(layerIdx).matID = 4;
    physLayers(layerIdx).nz = max(1, round(10 / opts.target_dz));
    z_current = z_current + 10;

    % 400 nm In0.08Ga0.92N
    layerIdx = layerIdx + 1;
    physLayers(layerIdx).name = 'InGaN_400nm';
    physLayers(layerIdx).z0 = z_current;
    physLayers(layerIdx).thickness = 400;
    physLayers(layerIdx).matID = 5;
    physLayers(layerIdx).nz = max(2, round(400 / opts.target_dz));
    z_current = z_current + 400;

    % 30 nm In0.08Ga0.92N
    layerIdx = layerIdx + 1;
    physLayers(layerIdx).name = 'InGaN_30nm_1';
    physLayers(layerIdx).z0 = z_current;
    physLayers(layerIdx).thickness = 30;
    physLayers(layerIdx).matID = 6;
    physLayers(layerIdx).nz = max(1, round(30 / opts.target_dz));
    z_current = z_current + 30;

    % 2.6 nm In0.18Ga0.82N
    layerIdx = layerIdx + 1;
    physLayers(layerIdx).name = 'InGaN_2p6nm';
    physLayers(layerIdx).z0 = z_current;
    physLayers(layerIdx).thickness = 2.6;
    physLayers(layerIdx).matID = 7;
    physLayers(layerIdx).nz = 1;
    z_current = z_current + 2.6;

    % 30 nm In0.08Ga0.92N (top)
    layerIdx = layerIdx + 1;
    physLayers(layerIdx).name = 'InGaN_30nm_2';
    physLayers(layerIdx).z0 = z_current;
    physLayers(layerIdx).thickness = 30;
    physLayers(layerIdx).matID = 8;
    physLayers(layerIdx).nz = max(1, round(30 / opts.target_dz));
    z_current = z_current + 30;

    zTop = z_current;

    % Clamp zBank if needed (should not exceed pillar top)
    if opts.zBank > zTop
        warning('zBank %.1f exceeds pillar top %.1f, clamping to zTop', opts.zBank, zTop);
        opts.zBank = zTop;
    end

    % Now insert interface sublayers between all adjacent physical layers
    % AND ensure z=0 and zBank are slab boundaries
    % Split layers that straddle critical boundaries
    layers = struct([]);
    layerCount = 0;

    criticalBoundaries = [0, opts.zBank];

    for i = 1:length(physLayers)
        z_bot = physLayers(i).z0;
        z_top = physLayers(i).z0 + physLayers(i).thickness;

        % Find critical boundaries within this layer
        boundaries_in_layer = criticalBoundaries(criticalBoundaries > z_bot & criticalBoundaries < z_top);

        % Build list of split points
        splitPoints = [z_bot, boundaries_in_layer, z_top];

        % Create sub-layers for each segment
        for j = 1:length(splitPoints)-1
            layerCount = layerCount + 1;
            z_start = splitPoints(j);
            z_end = splitPoints(j+1);
            thickness_segment = z_end - z_start;

            % Calculate nz for this segment (proportional to thickness)
            nz_segment = max(1, round(physLayers(i).nz * thickness_segment / physLayers(i).thickness));

            if length(splitPoints) > 2
                % This is a split layer
                layers(layerCount).name = sprintf('%s_seg%d', physLayers(i).name, j);
            else
                % Not split
                layers(layerCount).name = physLayers(i).name;
            end

            layers(layerCount).z0 = z_start;
            layers(layerCount).thickness = thickness_segment;
            layers(layerCount).matID = physLayers(i).matID;
            layers(layerCount).nz = nz_segment;
            layers(layerCount).isInterface = false;
        end

        % Add interface after this layer (if not last)
        if i < length(physLayers)
            layerCount = layerCount + 1;
            layers(layerCount).name = sprintf('Interface_%d_%d', i, i+1);
            layers(layerCount).z0 = z_top;
            layers(layerCount).thickness = opts.tInterface;
            layers(layerCount).matID = 100 + i;
            layers(layerCount).nz = 1;
            layers(layerCount).isInterface = true;
        end
    end

    % Build z-levels array from layers
    zLevels = layers(1).z0;
    for i = 1:length(layers)
        % Add internal divisions for this layer
        z_start = layers(i).z0;
        z_end = z_start + layers(i).thickness;
        nz = layers(i).nz;

        % Create internal z-levels
        z_layer = linspace(z_start, z_end, nz + 1);
        zLevels = [zLevels, z_layer(2:end)];
    end

    % Remove duplicates and sort
    zLevels = unique(zLevels);
    zLevels = sort(zLevels);

    % Force exact z=0 and zBank
    [~, idx0] = min(abs(zLevels));
    if abs(zLevels(idx0)) < 1e-6
        zLevels(idx0) = 0;
    end

    [~, idxB] = min(abs(zLevels - opts.zBank));
    if abs(zLevels(idxB) - opts.zBank) < 1e-6
        zLevels(idxB) = opts.zBank;
    end
end

function sf = createShapeFunction(sfName)
    % Create appropriate shape function object
    switch upper(sfName)
        case 'H27'
            sf = ShapeFunctionH27();
        case 'H8'
            sf = ShapeFunctionH8();
        otherwise
            error('Unknown shape function: %s. Use H27 or H8.', sfName);
    end
end

function r_pillar = pillarRadius(z, rTop, zTop, alpha_rad)
    % Compute pillar radius at height z (cone expanding downward)
    % r_pillar(z) = rTop + (zTop - z) * tan(alpha)
    r_pillar = rTop + (zTop - z) * tan(alpha_rad);
end

function r_boundary = computeBoundary(z, opts, zTop)
    % Compute radial boundary r_boundary(z) for solid/void interface
    % This represents the etched surface profile
    %
    % For z > 0:
    %   - Inside pillar cone: r = pillar radius
    %   - Outside pillar: depends on bank height
    %
    % For z <= 0:
    %   - Smooth trench profile with bank

    alpha_rad = opts.alpha_deg * pi / 180;
    r_pillar_z = pillarRadius(z, opts.rTop, zTop, alpha_rad);

    if z > opts.zBank
        % Above bank: only pillar exists
        r_boundary = r_pillar_z;
    elseif z > 0
        % Between z=0 and zBank: pillar + bank region
        % Bank extends to some radius beyond pillar
        % Use smooth transition
        r_bank_max = r_pillar_z + opts.trenchWidth * 0.5;
        r_boundary = r_bank_max;
    else
        % z <= 0: trench region
        % Smooth profile: starts at some radius at z=0, deepens, then returns to pillar base
        % Use a smooth function (cosine-based)

        % At z=0, boundary is near pillar base radius
        r_pillar_0 = pillarRadius(0, opts.rTop, zTop, alpha_rad);
        r_trench_center = r_pillar_0 + opts.trenchWidth * 0.5;

        % Depth profile: smooth dip
        depth_factor = smoothstep(max(0, min(1, -z / opts.trenchDepth)));
        r_trench = r_trench_center + opts.trenchWidth * 0.3 * depth_factor;

        % At very bottom, return to pillar radius
        if z < -opts.trenchDepth
            blend = smoothstep(max(0, min(1, (-z - opts.trenchDepth) / opts.Hcurve)));
            r_boundary = (1 - blend) * r_trench + blend * r_pillar_z;
        else
            r_boundary = r_trench;
        end
    end
end

function s = smoothstep(x)
    % Smooth interpolation function (3x^2 - 2x^3)
    s = x .* x .* (3 - 2 * x);
end

function isInside = isPointInSolid(x, y, z, opts, zTop)
    % Determine if point (x,y,z) is inside solid material
    % (not in etched/void region)

    r = sqrt(x^2 + y^2);
    alpha_rad = opts.alpha_deg * pi / 180;
    r_pillar_z = pillarRadius(z, opts.rTop, zTop, alpha_rad);

    if z > opts.zBank
        % Above bank: only inside pillar cone
        isInside = (r <= r_pillar_z);
    elseif z > 0
        % Between z=0 and zBank
        % Inside if within bank radius
        r_boundary_z = computeBoundary(z, opts, zTop);
        isInside = (r <= r_boundary_z);
    else
        % z <= 0: check against trench boundary
        r_boundary_z = computeBoundary(z, opts, zTop);
        isInside = (r <= r_boundary_z);
    end
end

function block = generateCylindricalBlock(Rin, Rout, zLevels, nr, ntheta, ...
                                          layers, opts, sf, zTop)
    % Generate cylindrical annular block
    % Quarter circle: theta from 0 to pi/2

    % Radial divisions
    r_edges = linspace(Rin, Rout, nr + 1);

    % Angular divisions (quarter circle)
    theta_edges = linspace(0, pi/2, ntheta + 1);

    % Generate nodes
    nodes = [];
    nodeMap = zeros(length(r_edges), length(theta_edges), length(zLevels));
    nodeIdx = 1;

    for iz = 1:length(zLevels)
        z = zLevels(iz);
        for it = 1:length(theta_edges)
            theta = theta_edges(it);
            for ir = 1:length(r_edges)
                r = r_edges(ir);
                x = r * cos(theta);
                y = r * sin(theta);

                nodes(nodeIdx, :) = [x, y, z];
                nodeMap(ir, it, iz) = nodeIdx;
                nodeIdx = nodeIdx + 1;
            end
        end
    end

    % Generate elements
    elems = [];
    matID = [];

    for iz = 1:length(zLevels)-1
        % Find which layer this z-slab belongs to
        z_mid = (zLevels(iz) + zLevels(iz+1)) / 2;
        layerIdx = findLayer(z_mid, layers);

        for it = 1:ntheta
            for ir = 1:nr
                % Get 8 corner nodes for this hex
                corners = [
                    nodeMap(ir,   it,   iz)
                    nodeMap(ir+1, it,   iz)
                    nodeMap(ir+1, it+1, iz)
                    nodeMap(ir,   it+1, iz)
                    nodeMap(ir,   it,   iz+1)
                    nodeMap(ir+1, it,   iz+1)
                    nodeMap(ir+1, it+1, iz+1)
                    nodeMap(ir,   it+1, iz+1)
                ];

                % Check if element is in solid region
                % Sample at element centroid
                x_center = mean(nodes(corners, 1));
                y_center = mean(nodes(corners, 2));
                z_center = mean(nodes(corners, 3));

                inSolid = isPointInSolid(x_center, y_center, z_center, opts, zTop);

                if inSolid || strcmp(opts.etch_mode, 'air_material')
                    % Build full element connectivity using shape function
                    elemNodes = buildElementConnectivity(corners, nodeMap, ...
                        ir, it, iz, sf);

                    elems(end+1, :) = elemNodes;

                    if inSolid
                        matID(end+1) = layers(layerIdx).matID;
                    else
                        matID(end+1) = 999; % Air material
                    end
                end
            end
        end
    end

    % Create mesh block structure
    block.nodes = nodes;
    block.elems = elems;
    block.matID = matID(:);
end

function block = generateTransitionBlock(Rin, Rout, zLevels, nr, ntheta, ...
                                         layers, opts, sf, zTop)
    % Generate transition block with circle-to-square warping
    % Inner boundary at r=Rin (circular)
    % Outer boundary at square [0, Rout] x [0, Rout]

    % Radial parameter s from 0 (inner circle) to 1 (outer square)
    s_edges = linspace(0, 1, nr + 1);

    % Angular divisions
    theta_edges = linspace(0, pi/2, ntheta + 1);

    % Generate nodes with warping
    nodes = [];
    nodeMap = zeros(length(s_edges), length(theta_edges), length(zLevels));
    nodeIdx = 1;

    for iz = 1:length(zLevels)
        z = zLevels(iz);
        for it = 1:length(theta_edges)
            theta = theta_edges(it);
            for is = 1:length(s_edges)
                s = s_edges(is);

                % Warp from circle to square
                if s == 0
                    % Inner circle
                    r = Rin;
                    x = r * cos(theta);
                    y = r * sin(theta);
                else
                    % Blend between circle and square
                    r_circ = Rin + s * (Rout - Rin);
                    x_circ = r_circ * cos(theta);
                    y_circ = r_circ * sin(theta);

                    % Square corner
                    if theta <= pi/4
                        x_sq = Rout;
                        y_sq = Rout * tan(theta);
                    else
                        x_sq = Rout / tan(theta);
                        y_sq = Rout;
                    end

                    % Blend with smoothstep
                    blend = smoothstep(s);
                    x = (1 - blend) * x_circ + blend * x_sq;
                    y = (1 - blend) * y_circ + blend * y_sq;
                end

                nodes(nodeIdx, :) = [x, y, z];
                nodeMap(is, it, iz) = nodeIdx;
                nodeIdx = nodeIdx + 1;
            end
        end
    end

    % Generate elements (similar to cylindrical)
    elems = [];
    matID = [];

    for iz = 1:length(zLevels)-1
        z_mid = (zLevels(iz) + zLevels(iz+1)) / 2;
        layerIdx = findLayer(z_mid, layers);

        for it = 1:ntheta
            for is = 1:nr
                corners = [
                    nodeMap(is,   it,   iz)
                    nodeMap(is+1, it,   iz)
                    nodeMap(is+1, it+1, iz)
                    nodeMap(is,   it+1, iz)
                    nodeMap(is,   it,   iz+1)
                    nodeMap(is+1, it,   iz+1)
                    nodeMap(is+1, it+1, iz+1)
                    nodeMap(is,   it+1, iz+1)
                ];

                x_center = mean(nodes(corners, 1));
                y_center = mean(nodes(corners, 2));
                z_center = mean(nodes(corners, 3));

                inSolid = isPointInSolid(x_center, y_center, z_center, opts, zTop);

                if inSolid || strcmp(opts.etch_mode, 'air_material')
                    elemNodes = buildElementConnectivity(corners, nodeMap, ...
                        is, it, iz, sf);

                    elems(end+1, :) = elemNodes;

                    if inSolid
                        matID(end+1) = layers(layerIdx).matID;
                    else
                        matID(end+1) = 999;
                    end
                end
            end
        end
    end

    block.nodes = nodes;
    block.elems = elems;
    block.matID = matID(:);
end

function elemNodes = buildElementConnectivity(corners, nodeMap, i, j, k, sf)
    % Build element connectivity from corner nodes using shape function
    % corners: [8 x 1] corner node indices
    % nodeMap: 3D array of node indices
    % i, j, k: element position in structured grid
    % sf: ShapeFunction object

    % Get local node positions from shape function
    localNodes = sf.localNodes; % [nen x 3] in local coords [-1, 1]
    nen = size(localNodes, 1);
    elemNodes = zeros(1, nen);

    % Map local nodes to global connectivity
    % Local coords: xi, eta, zeta in [-1, 1]
    % Map to structured indices

    for n = 1:nen
        xi = localNodes(n, 1);
        eta = localNodes(n, 2);
        zeta = localNodes(n, 3);

        % Map to structured grid indices
        % xi:   -1 -> i,   +1 -> i+1
        % eta:  -1 -> j,   +1 -> j+1
        % zeta: -1 -> k,   +1 -> k+1

        ii = i + (1 + xi) / 2;
        jj = j + (1 + eta) / 2;
        kk = k + (1 + zeta) / 2;

        % Round to nearest integer (for edge/face nodes)
        ii = round(ii);
        jj = round(jj);
        kk = round(kk);

        elemNodes(n) = nodeMap(ii, jj, kk);
    end
end

function layerIdx = findLayer(z, layers)
    % Find which layer contains height z
    for i = 1:length(layers)
        z0 = layers(i).z0;
        z1 = z0 + layers(i).thickness;
        if z >= z0 && z <= z1
            layerIdx = i;
            return;
        end
    end
    % If not found, return last layer (shouldn't happen)
    layerIdx = length(layers);
end

function verifyMesh(mesh, opts, zTop, zLevels, layers)
    % Verify mesh properties
    fprintf('=== Mesh Verification ===\n');

    % Check interface layers
    nInterfaces = sum([layers.isInterface]);
    fprintf('Interface sublayers: %d\n', nInterfaces);

    for i = 1:length(layers)
        if layers(i).isInterface
            assert(layers(i).nz == 1, 'Interface layer %s has nz=%d, expected 1', ...
                layers(i).name, layers(i).nz);
        end
    end
    fprintf('All interface layers have nz=1: PASS\n');

    % Check maximum z outside pillar
    alpha_rad = opts.alpha_deg * pi / 180;

    if ~isempty(mesh.nodes)
        maxZ_outside = -inf;
        for i = 1:size(mesh.nodes, 1)
            x = mesh.nodes(i, 1);
            y = mesh.nodes(i, 2);
            z = mesh.nodes(i, 3);
            r = sqrt(x^2 + y^2);

            if z > 0
                r_pillar_z = pillarRadius(z, opts.rTop, zTop, alpha_rad);
                if r > r_pillar_z + 1e-6 % Outside pillar
                    maxZ_outside = max(maxZ_outside, z);
                end
            end
        end

        if maxZ_outside > -inf
            fprintf('Max z outside pillar: %.2f nm (expected ~%.2f nm)\n', ...
                maxZ_outside, opts.zBank);
        else
            fprintf('No nodes found outside pillar for z>0\n');
        end
    end

    % Check node z-coordinates are from zLevels (no warping)
    if ~isempty(mesh.nodes)
        uniqueZ = unique(mesh.nodes(:, 3));
        fprintf('Unique z-coordinates in mesh: %d\n', length(uniqueZ));
        fprintf('z-levels in stack: %d\n', length(zLevels));

        % All mesh z should match zLevels
        allMatch = true;
        for i = 1:length(uniqueZ)
            [minDist, ~] = min(abs(zLevels - uniqueZ(i)));
            if minDist > 1e-6
                allMatch = false;
                fprintf('WARNING: z=%.4f does not match any zLevel\n', uniqueZ(i));
            end
        end
        if allMatch
            fprintf('All node z-coordinates match zLevels: PASS\n');
        end
    end

    fprintf('=========================\n\n');
end

% EXAMPLEFEMUSAGE - Demonstrate how to use pillar mesh in FEM analysis
%
% This script shows a complete workflow from mesh generation to
% finite element model setup (but does not solve, just demonstrates setup)

clear; clc;
fprintf('\n=== PILLAR MESH FEM USAGE EXAMPLE ===\n\n');

%% Step 1: Generate the mesh
fprintf('Step 1: Generating mesh...\n');

opts = struct();
opts.target_dz = 20;       % 20 nm vertical resolution
opts.nr_core = 8;          % 8 radial divisions in core
opts.nr_trench = 8;        % 8 radial divisions in trench
opts.nr_outer = 6;         % 6 radial divisions in outer zone
opts.ntheta = 8;           % 8 angular divisions
opts.sfName = 'H27';       % 27-node hexahedral elements

mesh = generatePillarQuarterMesh(opts);

fprintf('  Generated mesh with %d nodes and %d elements\n', ...
    size(mesh.nodes, 1), size(mesh.elems, 1));
fprintf('\n');

%% Step 2: Create finite element object
fprintf('Step 2: Creating finite element object...\n');

sf = mesh.sf;  % Use the shape function from the mesh
fe = SolidElasticElem(sf, mesh.elems);

fprintf('  Created FE object with shape function: %s\n', class(sf));
fprintf('\n');

%% Step 3: Define materials
fprintf('Step 3: Defining materials...\n');

% Material properties (example values - adjust for your materials)
materials = struct();

% Material 1: Extra substrate (if used)
materials(1).name = 'Substrate_extra';
materials(1).E = 300e3;    % Young's modulus [MPa]
materials(1).nu = 0.25;    % Poisson's ratio

% Material 2: GaN substrate
materials(2).name = 'GaN_substrate';
materials(2).E = 295e3;    % Young's modulus [MPa]
materials(2).nu = 0.183;   % Poisson's ratio

% Material 3: Porous GaN:Si
materials(3).name = 'GaN_porous';
materials(3).E = 150e3;    % Reduced E due to porosity
materials(3).nu = 0.20;

% Material 4: GaN (10 nm pillar layer)
materials(4).name = 'GaN_10nm';
materials(4).E = 295e3;
materials(4).nu = 0.183;

% Material 5: In0.08Ga0.92N (400 nm)
materials(5).name = 'InGaN_400nm';
materials(5).E = 280e3;    % Approximate value
materials(5).nu = 0.20;

% Material 6: In0.08Ga0.92N (30 nm, first)
materials(6).name = 'InGaN_30nm_1';
materials(6).E = 280e3;
materials(6).nu = 0.20;

% Material 7: In0.18Ga0.82N (2.6 nm quantum well)
materials(7).name = 'InGaN_2p6nm_QW';
materials(7).E = 250e3;    % Lower E with more In
materials(7).nu = 0.22;

% Material 8: In0.08Ga0.92N (30 nm, top)
materials(8).name = 'InGaN_30nm_2';
materials(8).E = 280e3;
materials(8).nu = 0.20;

% Display material assignments
uniqueMats = unique(mesh.matID);
fprintf('  Material assignments:\n');
for i = 1:length(uniqueMats)
    matID = uniqueMats(i);
    nElems = sum(mesh.matID == matID);

    if matID <= length(materials)
        fprintf('    Mat %d (%s): %d elements\n', ...
            matID, materials(matID).name, nElems);
    elseif matID >= 100 && matID < 200
        fprintf('    Mat %d (Interface_%d): %d elements\n', ...
            matID, matID-100, nElems);
    elseif matID == 999
        fprintf('    Mat %d (Air): %d elements\n', matID, nElems);
    end
end
fprintf('\n');

%% Step 4: Assign materials to elements
fprintf('Step 4: Assigning materials to finite elements...\n');

% For this example, we'll create a single "averaged" material
% In practice, you would assign different materials to different element groups

% Example: Create a single material for all solid elements
avgMaterial = SolidMaterial('pillar_avg');
avgMaterial.setElasticIzo(280e3, 0.20);  % Average InGaN properties
avgMaterial.setElasticIzoGrad();         % Set up gradient

% Assign to solid elements (matID != 999)
solidElemIndices = find(mesh.matID ~= 999);
fe.setMaterial(avgMaterial);

fprintf('  Assigned averaged material to %d solid elements\n', length(solidElemIndices));
fprintf('  (In practice, assign different materials based on mesh.matID)\n');
fprintf('\n');

%% Step 5: Setup analysis (example - no solve)
fprintf('Step 5: Setting up analysis framework...\n');

% Create analysis object
analysis = LinearElasticityWeighted(fe, mesh, false);

% Define boundary conditions (example: fix bottom, load top)

% Fix bottom surface (z = zmin)
zmin = min(mesh.nodes(:,3));
fixedSelector = Selector(@(x)( abs(x(:,3) - zmin) < 1e-6 ));
analysis.fixNodes(fixedSelector, ["ux" "uy" "uz"]);

% Apply load to pillar top (example: uniform pressure)
zmax = max(mesh.nodes(:,3));
loadSelector = Selector(@(x)( abs(x(:,3) - zmax) < 1e-6 ));
pressure = 10;  % MPa

% Note: Surface load requires identifying top surface elements
% This is simplified - actual implementation would use element faces
% analysis.elementLoadSurfaceIntegral("global", loadSelector, ...
%     ["ux" "uy" "uz"], @(x)( x*0 + [0 0 -pressure] ));

fprintf('  Boundary conditions:\n');
fprintf('    Fixed nodes at z = %.2f nm\n', zmin);
fprintf('    Load region at z = %.2f nm (pressure = %.1f MPa)\n', zmax, pressure);
fprintf('\n');

%% Step 6: Print problem info
fprintf('Step 6: Problem information...\n');
analysis.printProblemInfo();
fprintf('\n');

%% Step 7: Visualization
fprintf('Step 7: Visualizing mesh and materials...\n');

% Visualize the mesh
visualizePillarMesh(mesh);

% Add slice view at quantum well location
% The 2.6 nm quantum well is approximately at z ≈ 440 nm
% (10 + 400 + 30 = 440 nm above z=0)

figure('Name', 'Quantum Well Region', 'Position', [100, 100, 800, 600]);
qw_z = 440;  % Approximate QW location
tol = 5;
idx_qw = abs(mesh.nodes(:,3) - qw_z) < tol;

if sum(idx_qw) > 0
    scatter(mesh.nodes(idx_qw, 1), mesh.nodes(idx_qw, 2), 20, 'r', 'filled');
    axis equal; grid on;
    xlabel('x [nm]'); ylabel('y [nm]');
    title(sprintf('Nodes near Quantum Well (z ≈ %.0f ± %.0f nm)', qw_z, tol));
end

fprintf('  Created visualization figures\n');
fprintf('\n');

%% Summary
fprintf('=== SETUP COMPLETE ===\n');
fprintf('The mesh is ready for FEM analysis.\n');
fprintf('To solve:\n');
fprintf('  1. Refine material assignments using mesh.matID\n');
fprintf('  2. Apply appropriate boundary conditions for your problem\n');
fprintf('  3. Call analysis.solveWeighted() or similar\n');
fprintf('  4. Post-process results with analysis.plotMaps(), etc.\n');
fprintf('\n');
fprintf('Key mesh properties:\n');
fprintf('  mesh.nodes  : %d x 3 array of node coordinates\n', size(mesh.nodes, 1));
fprintf('  mesh.elems  : %d x %d element connectivity\n', size(mesh.elems, 1), size(mesh.elems, 2));
fprintf('  mesh.matID  : %d x 1 material IDs\n', length(mesh.matID));
fprintf('  mesh.sf     : %s shape function\n', class(mesh.sf));
fprintf('\n');

%% Helper: Export material map
fprintf('Exporting material map...\n');

% Create element centroid coordinates with material IDs
nElems = size(mesh.elems, 1);
elemData = zeros(nElems, 4);  % [x, y, z, matID]

for i = 1:nElems
    elemData(i, 1:3) = mean(mesh.nodes(mesh.elems(i, :), :), 1);
    elemData(i, 4) = mesh.matID(i);
end

% Save to file (optional)
% save('pillar_element_data.mat', 'elemData');
fprintf('Element data prepared (not saved)\n');
fprintf('  Format: [x_centroid, y_centroid, z_centroid, materialID]\n');
fprintf('\n');

fprintf('=== EXAMPLE COMPLETE ===\n');

% TESTPILLARMESH - Demo and test script for pillar mesh generator
%
% This script demonstrates usage of generatePillarQuarterMesh and
% verifies the mesh properties.

%% Setup
clear; clc;
fprintf('\n=== PILLAR MESH TEST ===\n\n');

%% Test 1: Default configuration
fprintf('TEST 1: Default configuration\n');
fprintf('-----------------------------\n');

opts = struct();
opts.target_dz = 25; % Slightly coarser for faster demo

mesh1 = generatePillarQuarterMesh(opts);

fprintf('\n');

%% Test 2: Higher resolution
fprintf('TEST 2: Higher resolution\n');
fprintf('-------------------------\n');

opts2 = struct();
opts2.target_dz = 15;
opts2.nr_core = 12;
opts2.nr_trench = 10;
opts2.nr_outer = 8;
opts2.ntheta = 12;

% For demo, skip to save time (uncomment to run)
% mesh2 = generatePillarQuarterMesh(opts2);

fprintf('Skipped for demo (uncomment to run)\n\n');

%% Test 3: Using H8 elements
fprintf('TEST 3: H8 elements (linear hexahedra)\n');
fprintf('--------------------------------------\n');

opts3 = struct();
opts3.sfName = 'H8';
opts3.target_dz = 30;

mesh3 = generatePillarQuarterMesh(opts3);

fprintf('\n');

%% Test 4: Air material mode (optional)
fprintf('TEST 4: Air material mode\n');
fprintf('-------------------------\n');

opts4 = struct();
opts4.etch_mode = 'air_material';
opts4.target_dz = 30;

% Uncomment to run
% mesh4 = generatePillarQuarterMesh(opts4);

fprintf('Skipped for demo (uncomment to run)\n\n');

%% Visualize mesh
fprintf('Visualization\n');
fprintf('-------------\n');

if size(mesh1.nodes, 1) < 100000
    % Use custom visualization function
    visualizePillarMesh(mesh1);
    fprintf('Created comprehensive visualization\n');
else
    fprintf('Mesh too large for visualization (%d nodes)\n', size(mesh1.nodes, 1));
end

%% Print summary statistics
fprintf('\n');
fprintf('=== SUMMARY STATISTICS ===\n');
fprintf('Mesh 1 (H27, default):\n');
fprintf('  Nodes: %d\n', size(mesh1.nodes, 1));
fprintf('  Elements: %d\n', size(mesh1.elems, 1));
fprintf('  Element type: %s\n', class(mesh1.sf));
fprintf('  Nodes per element: %d\n', size(mesh1.sf.localNodes, 1));
fprintf('\n');

fprintf('Mesh 3 (H8):\n');
fprintf('  Nodes: %d\n', size(mesh3.nodes, 1));
fprintf('  Elements: %d\n', size(mesh3.elems, 1));
fprintf('  Element type: %s\n', class(mesh3.sf));
fprintf('  Nodes per element: %d\n', size(mesh3.sf.localNodes, 1));
fprintf('\n');

%% Verify mesh quality
fprintf('=== MESH QUALITY CHECKS ===\n');

% Check 1: All element nodes are valid
maxNodeIdx = max(mesh1.elems(:));
if maxNodeIdx <= size(mesh1.nodes, 1)
    fprintf('All element nodes valid: PASS\n');
else
    fprintf('Invalid element nodes found: FAIL\n');
end

% Check 2: No duplicate nodes in elements
hasDuplicates = false;
for i = 1:min(100, size(mesh1.elems, 1)) % Check first 100 elements
    if length(unique(mesh1.elems(i, :))) < length(mesh1.elems(i, :))
        hasDuplicates = true;
        break;
    end
end
if ~hasDuplicates
    fprintf('No duplicate nodes in elements: PASS\n');
else
    fprintf('Duplicate nodes found in elements: FAIL\n');
end

% Check 3: Material IDs assigned
if ~isempty(mesh1.matID) && length(mesh1.matID) == size(mesh1.elems, 1)
    fprintf('Material IDs assigned: PASS\n');
    uniqueMats = unique(mesh1.matID);
    fprintf('  Unique material IDs: %s\n', mat2str(uniqueMats'));
else
    fprintf('Material ID issue: FAIL\n');
end

% Check 4: Node coordinate ranges
fprintf('\nNode coordinate ranges:\n');
fprintf('  x: [%.2f, %.2f] nm\n', min(mesh1.nodes(:,1)), max(mesh1.nodes(:,1)));
fprintf('  y: [%.2f, %.2f] nm\n', min(mesh1.nodes(:,2)), max(mesh1.nodes(:,2)));
fprintf('  z: [%.2f, %.2f] nm\n', min(mesh1.nodes(:,3)), max(mesh1.nodes(:,3)));

fprintf('\n=== TEST COMPLETE ===\n');

figure;
sf = ShapeFunctionH27();
fe = SolidElasticElem( sf, mesh1.elems );
fe.plot(mesh1.nodes);
       

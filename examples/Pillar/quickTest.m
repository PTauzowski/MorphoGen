% Quick test to verify mesh generator runs
clear; clc;

fprintf('Quick test of pillar mesh generator...\n');

% Very coarse mesh for quick testing
opts = struct();
opts.target_dz = 50;  % Coarse
opts.nr_core = 4;
opts.nr_trench = 4;
opts.nr_outer = 3;
opts.ntheta = 4;
opts.sfName = 'H8';  % Use H8 for speed

try
    mesh = generatePillarQuarterMesh(opts);
    fprintf('SUCCESS! Generated mesh with %d nodes, %d elements\n', ...
        size(mesh.nodes, 1), size(mesh.elems, 1));
catch ME
    fprintf('ERROR: %s\n', ME.message);
    fprintf('Stack:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

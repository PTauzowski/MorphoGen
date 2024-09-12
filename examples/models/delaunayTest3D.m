% Create a PDE model
model = createpde;

% Define the geometry: a unit cube with dimensions 10x10x10
gm = multicuboid(10, 10, 10);

% Assign the geometry to the PDE model
model.Geometry = gm;

% Define the mesh size function
meshSizeFcn = @(location, state) 1 ./ (1 + sin(location.x) .* cos(location.y) .* exp(-location.z / 5));

% Generate the mesh with custom density function
generateMesh(model, 'Hmax', meshSizeFcn);

% Extract the mesh
mesh = model.Mesh;

% Plot the mesh
figure;
pdemesh(mesh);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Smooth Tetrahedral Mesh with Density Function');


clear;
close all;

Ecl  = 25.0E9; 
Ebm  = 25.0E9;
Egr  = 80.0E9;
nugr = 0.3;

b=0.25;
hbm=0.8;
hcl=0.7;
lspan=5;
hfloor=3;
nspan=3;
nfloor=3;


lgr=30;
hgr=15;
resgr=5;

%model = FrameOnElasticGroundModel( Ecl, hcl, Ebm, hbm, b, lspan, hfloor, nspan, nfloor, Egr, nugr, lgr, hgr, resgr, ShapeFunctionL4 );


xp=[(lgr-nspan*lspan)/2 hgr];

shapeFn2D=ShapeFunctionL16();
       
mesh=Mesh();
frameElem = Frame2D( mesh.addHframe(nspan,lspan,nfloor,hfloor,xp));
planeElem = PlaneStressElem( shapeFn2D, mesh.addRectMeshArray2D(0,0,[xp(1) repelem(lspan,1,nspan) xp(1)], xp(2), resgr, shapeFn2D.pattern) );

model = FEModel( { frameElem planeElem }, mesh );
model.plotWired();
edgeX0 = selectX(0.0);
model.fixDOF( edgeX0, ["ux" "uy"] , [0 0]);
%model.plotSelectedNodes(edgeX0,'m');
model.plotNodes(".",'b');


% Parameters
N = 10;  % Number of intervals along each axis
L = 1;   % Length of the cube along each axis

% Generate grid points
[x, y, z] = meshgrid(linspace(0, L, N + 1), linspace(0, L, N + 1), linspace(0, L, N + 1));

% Reshape the grid points for each axis
% X-axis edges
X1 = reshape([x(:, 1:end-1, :); x(:, 2:end, :)], 1, []);  % Repeat x along X edges
Y1 = reshape([y(:, 1:end-1, :); y(:, 2:end, :)], 1, []);
Z1 = reshape([z(:, 1:end-1, :); z(:, 2:end, :)], 1, []);

% Y-axis edges
X2 = reshape([x(1:end-1, :, :); x(2:end, :, :)], 1, []);  % Repeat y along Y edges
Y2 = reshape([y(1:end-1, :, :); y(2:end, :, :)], 1, []);
Z2 = reshape([z(1:end-1, :, :); z(2:end, :, :)], 1, []);

% Z-axis edges
X3 = reshape([x(:, :, 1:end-1); x(:, :, 2:end)], 1, []);  % Repeat z along Z edges
Y3 = reshape([y(:, :, 1:end-1); y(:, :, 2:end)], 1, []);
Z3 = reshape([z(:, :, 1:end-1); z(:, :, 2:end)], 1, []);

% Concatenate all edges and add NaN separators between polylines
X = [X1 NaN X2 NaN X3 NaN];
Y = [Y1 NaN Y2 NaN Y3 NaN];
Z = [Z1 NaN Z2 NaN Z3 NaN];

% Plot the 3D cubic mesh using NaN-separated polylines
figure;
plot3(X, Y, Z, 'Color', 'b', 'LineWidth', 1);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['3D Cubic Mesh with ', num2str(N), ' Intervals']);
grid on;

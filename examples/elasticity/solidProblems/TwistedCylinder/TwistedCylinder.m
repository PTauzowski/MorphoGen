% TwistedCylinder.m
% Thin-wall annular cylinder clamped at the bottom, loaded by a pure twisting
% moment at the top.  Solves linear elasticity and visualises principal stress
% trajectories on the outer surface.
%
% Geometry and material are taken from armModelDefaults("thin") so the problem
% is geometrically consistent with the manipulator-arm study, but the model is
% built and solved here completely independently.

clear; close all; clc; clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..', '..', '..');
addpath(genpath(projectRoot));

resultRoot = fullfile(scriptDir, 'results');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

%% ---- Material and geometry -------------------------------------------------
arm    = armModelDefaults("thin");
E      = arm.E;         % 2.0e9 Pa
nu     = arm.nu;        % 0.35
R      = arm.R;         % 0.14 m  outer radius
r      = arm.r;         % 0.136 m inner radius
h      = 0.45; %arm.h_seg;     % 0.15 m  cylinder height
res_th = arm.res_th;    % 1       radial element layers

% Twisting moment: tangential force Pz applied at lever arm R
Ms = arm.Pz * R;        % [N·m]

% Mesh resolution — same formula as ManipulatorModel3D
resCirc = round(2*pi*R / (R-r) * res_th);          % circumferential
resLen  = max(1, round(h / (R-r) * res_th));        % axial
ShapeFn = ShapeFunctionL8();

fprintf('Cylinder: R=%.4f m, r=%.4f m, h=%.4f m, t=%.4f m\n', R, r, h, R-r);
fprintf('Mesh: resCirc=%d, resLen=%d, res_th=%d\n', resCirc, resLen, res_th);
fprintf('Ms = %.4g N*m\n', Ms);

%% ---- Mesh ------------------------------------------------------------------
% addManipulatorHalfSegment3D with alpha=0 generates a straight annular tube.
mesh = Mesh();
mesh.addManipulatorHalfSegment3D(r, R, h, 0, res_th, resCirc, resLen, ShapeFn.localNodes);
nNodes = size(mesh.nodes, 1);
nElems = size(mesh.elems, 1);
fprintf('Mesh: %d nodes, %d elements\n', nNodes, nElems);

%% ---- FE setup --------------------------------------------------------------
material = SolidMaterial('mat');
material.setElasticIzo(E, nu);
material.setElasticIzoGrad();

fe = SolidElasticElem(ShapeFn, mesh.elems);
fe.setMaterial(material);

analysis = LinearElasticityWeighted(fe, mesh, false);

%% ---- Boundary conditions ---------------------------------------------------
zmin = min(mesh.nodes(:,3));
zmax = max(mesh.nodes(:,3));
ztol = (zmax - zmin) * 1e-6;

bottomFlag = mesh.nodes(:,3) < zmin + ztol;
topFlag    = mesh.nodes(:,3) > zmax - ztol;
topIds     = find(topFlag);

% Clamp bottom face: fix all DOFs
analysis.fixNodes(Selector(bottomFlag), ["ux" "uy" "uz"]);
fprintf('Clamped %d bottom nodes (z = %.4g m)\n', nnz(bottomFlag), zmin);

% Torsion on top face: tangential traction t = Ms/J * [-(y-yc), (x-xc), 0]
J  = pi * (R^4 - r^4) / 2;   % polar moment of inertia of annular section
xc = 0;  yc = 0;              % centroid (by symmetry)

traction_fn = @(x) Ms/J * [-(x(:,2) - yc), (x(:,1) - xc), zeros(size(x,1),1)];
analysis.elementLoadSurfaceIntegral("global", Selector(topFlag), ["ux" "uy" "uz"], traction_fn);

% Scale nodal loads to exactly reproduce Ms (corrects surface-integration error)
Pn      = analysis.Pnodal(topIds, :);
xn      = mesh.nodes(topIds, :);
Ms_act  = sum(Pn(:,2).*(xn(:,1)-xc) - Pn(:,1).*(xn(:,2)-yc));
if abs(Ms_act) > eps
    sc = Ms / Ms_act;
    analysis.Pnodal(topIds, 1) = analysis.Pnodal(topIds, 1) * sc;
    analysis.Pnodal(topIds, 2) = analysis.Pnodal(topIds, 2) * sc;
end
fprintf('Torsion: Ms_target=%.4g N*m,  Ms_before_scale=%.4g N*m\n', Ms, Ms_act);

%% ---- Solve linear elasticity -----------------------------------------------
analysis.solveWeighted(ones(nElems, 1));
fprintf('Solved.  Max |u| = %.4e m\n', max(abs(analysis.qnodal(:))));

% Analytical twist angle for verification: phi = Ms * h / (G * J)
G       = E / (2*(1+nu));
phi_ref = Ms * h / (G * J);
% Top node displacement: u_x = -phi*y, u_y = +phi*x for rigid rotation
id_top  = topIds(1);
x_top   = mesh.nodes(id_top, 1);
y_top   = mesh.nodes(id_top, 2);
ux_ref  = -phi_ref * y_top;
uy_ref  =  phi_ref * x_top;
ux_fem  = analysis.qnodal(id_top, 1);
uy_fem  = analysis.qnodal(id_top, 2);
fprintf('Twist angle: phi_analytic=%.4e rad,  phi_FEM=%.4e rad\n', ...
    phi_ref, atan2(uy_fem, x_top));
fprintf('Top-node displacement: analytic=(%.4e, %.4e),  FEM=(%.4e, %.4e) m\n', ...
    ux_ref, uy_ref, ux_fem, uy_fem);

%% ---- Compute stresses ------------------------------------------------------
analysis.initializeResults();
analysis.computeElementResults();

% Nodal-averaged stress — column indices in fe.results.nodal.all.
%
% Important: SolidElasticElem's shear strain/stress B-matrix order is
% [yz, xz, xy], although results.names currently labels columns 10:12 as
% [sxy, syz, sxz].  Use the actual tensor order here; otherwise the
% cylindrical tau_theta_z field is assembled from swapped Cartesian shears
% and the principal trajectories become nonphysical.
S      = fe.results.nodal.all;
iSxx = 7; iSyy = 8; iSzz = 9; iSyz = 10; iSxz = 11; iSxy = 12; iHM = 13;

tau_ref = Ms * R / J;            % analytical shear stress at outer surface
sHM_ref = sqrt(3) * tau_ref;     % equivalent von Mises for pure shear
fprintf('Analytical tau_max = %.4g Pa  (sHM_equiv = %.4g Pa)\n', tau_ref, sHM_ref);
fprintf('FEM sHM_max        = %.4g Pa\n', max(S(:, iHM)));

%% ---- Principal stress directions on outer surface --------------------------
nodes = mesh.nodes;
rNode = sqrt(nodes(:,1).^2 + nodes(:,2).^2);

% Outer surface: nodes within 30% of wall thickness from R
outerTol  = (R - r) * 0.3;
outerMask = abs(rNode - R) < outerTol;
outerIds  = find(outerMask);

xO  = nodes(outerIds, 1);
yO  = nodes(outerIds, 2);
zO  = nodes(outerIds, 3);
thO = atan2(yO, xO);             % circumferential angle, in [-pi, pi]

sth = sin(thO);   cth = cos(thO);

% Cartesian stress components at outer nodes
sxx = S(outerIds, iSxx);   syy = S(outerIds, iSyy);   szz = S(outerIds, iSzz);
sxy = S(outerIds, iSxy);   syz = S(outerIds, iSyz);   sxz = S(outerIds, iSxz);

% Surface stress tensor in (e_theta, e_z) frame:
%   e_theta = [-sin(theta), cos(theta), 0]   (circumferential)
%   e_z     = [0, 0, 1]                      (axial)
%
%   S_tt = e_theta' * [S] * e_theta = sxx*sin^2(th) - 2*sxy*sin(th)*cos(th) + syy*cos^2(th)
%   S_zz_s = e_z' * [S] * e_z     = szz
%   S_tz = e_theta' * [S] * e_z   = -sxz*sin(th) + syz*cos(th)
S_tt   = sxx.*sth.^2 - 2*sxy.*sth.*cth + syy.*cth.^2;
S_zz_s = szz;
S_tz   = -sxz.*sth + syz.*cth;

% Principal angle in surface (from e_theta toward e_z)
% For pure torsion: S_tt=0, S_zz_s=0, S_tz=tau => alpha_p = pi/4 (45 deg)
alpha_p = 0.5 * atan2(2*S_tz, S_tt - S_zz_s);

fprintf('Mean principal angle on outer surface: %.2f deg  (expected 45 for pure torsion)\n', ...
    mean(alpha_p(~isnan(alpha_p))) * 180/pi);
fprintf('Std principal angle on outer surface:  %.2f deg  (near 0 away from end effects)\n', ...
    std(alpha_p(~isnan(alpha_p))) * 180/pi);
fprintf('Mean tau_theta_z / tau_ref on outer surface: %.4f\n', ...
    mean(S_tz(~isnan(S_tz))) / tau_ref);

%% ---- Stress trajectory integration ----------------------------------------
% Parameterisation: trajectory in (theta, z) on the outer cylinder.
%
% For a unit-arc-length step ds, the (theta, z) increments are:
%   Family 1 (principal direction 1): dtheta = cos(alpha_p)/R * ds, dz = sin(alpha_p) * ds
%   Family 2 (perpendicular in surface): dtheta = -sin(alpha_p)/R * ds, dz = cos(alpha_p) * ds
%
% This satisfies: R^2*(dtheta/ds)^2 + (dz/ds)^2 = 1  (unit arc length on cylinder).

% Build scattered interpolant for alpha_p, extended periodically in theta
% to avoid boundary artefacts when trajectories cross theta = +/-pi.
th_ext = [thO - 2*pi; thO; thO + 2*pi];
z_ext  = [zO;  zO;  zO];
a_ext  = [alpha_p; alpha_p; alpha_p];
F_alpha = scatteredInterpolant(th_ext, z_ext, a_ext, 'linear', 'linear');

% Step size: half of the smallest surface element dimension
ds      = 0.5 * min(2*pi*R/resCirc, h/resLen);
maxStep = ceil(2.0 * h / (ds * sin(pi/4)));   % enough to cross full height twice

% Seed: 24 equally spaced circumferential positions at z just above the clamp
nSeeds  = 24;
th0_vec = linspace(-pi, pi*(1 - 2/nSeeds), nSeeds)';
z0      = zmin + (zmax - zmin) * 0.02;   % 2% above clamped base

traj1 = cell(nSeeds, 1);
traj2 = cell(nSeeds, 1);
for k = 1:nSeeds
    traj1{k} = traceTrajectory(th0_vec(k), z0, 1, F_alpha, R, zmin, zmax, ds, maxStep);
    traj2{k} = traceTrajectory(th0_vec(k), z0, 2, F_alpha, R, zmin, zmax, ds, maxStep);
end

%% ---- Figure: principal stress trajectories ---------------------------------
fig1 = figure('Color', 'white', 'Units', 'normalized', 'Position', [0.05 0.05 0.55 0.78]);
hold on; axis off; daspect([1 1 1]);
view([30 22]);
light('Position', [-1 -2 5], 'Style', 'infinite');
light('Position', [ 1  1 2], 'Style', 'local');
lighting flat; material dull;

% Outer cylinder surface (analytical, smooth background)
[th_surf, z_surf] = meshgrid(linspace(-pi, pi, resCirc+1), linspace(zmin, zmax, resLen+1));
Xs = R * cos(th_surf);
Ys = R * sin(th_surf);
surf(Xs, Ys, z_surf, 'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.55);

% Inner cylinder surface
surf(r*cos(th_surf), r*sin(th_surf), z_surf, ...
    'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.40);

% Top and bottom annular caps
th_cap = linspace(-pi, pi, resCirc+1);
Xcap   = [R*cos(th_cap); r*cos(th_cap)];
Ycap   = [R*sin(th_cap); r*sin(th_cap)];
surf(Xcap, Ycap, zmin*ones(2,resCirc+1), ...
    'FaceColor', [0.70 0.70 0.70], 'EdgeColor', 'none', 'FaceAlpha', 0.65);
surf(Xcap, Ycap, zmax*ones(2,resCirc+1), ...
    'FaceColor', [0.70 0.70 0.70], 'EdgeColor', 'none', 'FaceAlpha', 0.65);

% Principal stress trajectories — family 1: tension principal (red)
h1 = [];
for k = 1:nSeeds
    pts = traj1{k};
    if size(pts,1) < 2, continue; end
    xv = R*cos(pts(:,1));   yv = R*sin(pts(:,1));   zv = pts(:,2);
    p = plot3(xv, yv, zv, 'r-', 'LineWidth', 1.8);
    if isempty(h1), h1 = p; end
end

% Principal stress trajectories — family 2: compression principal (blue)
h2 = [];
for k = 1:nSeeds
    pts = traj2{k};
    if size(pts,1) < 2, continue; end
    xv = R*cos(pts(:,1));   yv = R*sin(pts(:,1));   zv = pts(:,2);
    p = plot3(xv, yv, zv, 'b-', 'LineWidth', 1.8);
    if isempty(h2), h2 = p; end
end

if ~isempty(h1) && ~isempty(h2)
    legend([h1 h2], {'sigma_1 (tension)', 'sigma_2 (compression)'}, ...
        'Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
end

title({'Twisted Cylinder — Principal Stress Trajectories', ...
    sprintf('Ms = %.4g N*m,   R = %.3f m,   r = %.3f m,   h = %.3f m', Ms, R, r, h)}, ...
    'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(fig1, fullfile(resultRoot, 'trajectories.png'), 'Resolution', 200);
savefig(fig1, fullfile(resultRoot, 'trajectories.fig'));

%% ---- Figure: von Mises stress on outer surface -----------------------------
fig2 = figure('Color', 'white', 'Units', 'normalized', 'Position', [0.55 0.05 0.43 0.65]);
hold on; axis off; daspect([1 1 1]);
view([30 22]);
light('Position', [-1 -2 5], 'Style', 'infinite');
lighting flat; material dull;

% Interpolate HM stress onto smooth surface grid
F_sHM   = scatteredInterpolant(thO, zO, S(outerIds, iHM), 'linear', 'nearest');
sHM_grid = reshape(F_sHM(th_surf(:), z_surf(:)), size(th_surf));

surf(Xs, Ys, z_surf, sHM_grid, 'EdgeColor', 'none');
surf(Xcap, Ycap, zmin*ones(2,resCirc+1), ...
    'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
surf(Xcap, Ycap, zmax*ones(2,resCirc+1), ...
    'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');

colormap(jet); cb = colorbar;
cb.Label.String = 'Von Mises stress [Pa]';
title({'Twisted Cylinder — Von Mises Stress', sprintf('Ms = %.4g N*m', Ms)}, ...
    'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(fig2, fullfile(resultRoot, 'vonMises.png'), 'Resolution', 200);
savefig(fig2, fullfile(resultRoot, 'vonMises.fig'));

fprintf('\nAll results saved to: %s\n', resultRoot);

%% ---- Local functions -------------------------------------------------------

function pts = traceTrajectory(th0, z0, family, F_alpha, R, zmin, zmax, ds, maxStep)
% Integrate one principal stress trajectory on the outer cylindrical surface.
%
%   th0, z0  : seed point in (theta, z) parameterisation
%   family   : 1 = principal direction 1,  2 = perpendicular direction
%   F_alpha  : scattered interpolant for principal angle alpha_p(theta, z)
%   R        : outer cylinder radius (needed for dtheta/ds = cos/sin / R)
%   zmin,zmax: axial limits — stop when trajectory exits [zmin, zmax]
%   ds       : arc-length step size [m]
%   maxStep  : maximum number of integration steps
%
% Returns pts: [N x 2] array of (theta, z) pairs along the trajectory.

    pts  = zeros(maxStep + 1, 2);
    pts(1,:) = [th0, z0];
    n  = 1;
    th = th0;   z = z0;
    prev_d = [0, 1];   % initial direction hint: go upward in z

    for s = 1:maxStep
        a = F_alpha(th, z);
        if isnan(a), break; end

        if family == 1
            d = [cos(a) / R,  sin(a)];    % [dtheta/ds, dz/ds]
        else
            d = [-sin(a) / R,  cos(a)];
        end

        % Maintain direction consistency in the surface metric:
        % tangent = [R*dtheta/ds, dz/ds].
        if dot([R*d(1), d(2)], [R*prev_d(1), prev_d(2)]) < 0
            d = -d;
        end
        prev_d = d;

        th_new = th + ds * d(1);
        z_new  = z  + ds * d(2);

        if z_new < zmin || z_new > zmax
            break;
        end

        n = n + 1;
        th = th_new;
        z  = z_new;
        pts(n,:) = [th, z];
    end

    pts = pts(1:n, :);
end

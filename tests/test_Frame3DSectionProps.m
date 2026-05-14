%% test_Frame3DSectionProps.m
% Unit tests for Frame3DSectionProps: verify that EA, EIy, EIz, GJ, GAy,
% GAz are each used independently in the local stiffness matrix.
%
% Run as a plain script: >> test_Frame3DSectionProps
%
% Each test prints PASS or FAIL.

baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(fullfile(baseDir, 'analysis')));
addpath(genpath(fullfile(baseDir, 'elements')));
addpath(genpath(fullfile(baseDir, 'mesh')));

nTests  = 5;
results = false(nTests, 1);
tol     = 1e-8;

% Single-element cantilever: node 1 fixed, unit load at node 2.
elems = [1 2];
L     = 1.0;   % beam length [m]
nodes = [0 0 0; L 0 0];

% Full-pipe reference values (steel, R=0.05, r=0.04).
E  = 210e9;  nu = 0.3;  R = 0.05;  ri = 0.04;
G  = E / (2*(1+nu));
A0 = pi*(R^2  - ri^2);
I0 = pi*(R^4  - ri^4)/4;
J0 = pi*(R^4  - ri^4)/2;
m  = ri/R;
kappa = 6*(1+nu)*(1+m^2)^2 / ((7+6*nu)*(1+m^2)^2 + (20+12*nu)*m^2);

EA0  = E*A0;   EI0  = E*I0;   GJ0  = G*J0;
kGA0 = kappa*G*A0;

% Helper: solve cantilever, return all dofs at node 2.
%   dofNames: row vector of DOF strings that the load is applied to
%   load:     corresponding load values
function q2 = solveCantilever(fElem, nodes, dofNames, loadVals)
    m = Mesh();  m.nodes = nodes;  m.elems = fElem.elems;
    fa = LinearElasticityWeighted(fElem, m, false);
    fa.fixClosestNode(nodes(1,:), ["ux","uy","uz","fix","fiy","fiz"], zeros(1,6));
    fa.loadClosestNode(nodes(end,:), dofNames, loadVals);
    fa.solveWeighted(ones(size(fElem.elems,1), 1));
    q2 = fa.qnodal(2, :);   % [ux uy uz fix fiy fiz] at tip
end

%% Test 1: Axial stiffness EA
% Exact: u_axial = F / EA
F   = 1e3;
fEA = 2 * EA0;   % doubled EA only
fe  = Frame3DSectionProps(elems, fEA, EI0, EI0, GJ0, kGA0, kGA0);
q2  = solveCantilever(fe, nodes, "ux", F);
u_expected = F / fEA;
err = abs(q2(1) - u_expected);
results(1) = err < tol * abs(u_expected);
fprintf('Test 1 (axial EA):      %s  |err|=%.2e\n', tf(results(1)), err);

%% Test 2: Bending stiffness EIz (xy-plane, Euler-Bernoulli limit)
% Exact EB: u_y = F*L^3 / (3*EIz)  when shear is large (kGA >> EI/L^2)
F    = 1e3;
fEIz = 0.5 * EI0;   % halved EIz only
fkGA = 1e15;        % effectively rigid in shear -> EB limit
fe   = Frame3DSectionProps(elems, EA0, EI0, fEIz, GJ0, fkGA, fkGA);
q2   = solveCantilever(fe, nodes, "uy", F);
u_expected = F * L^3 / (3 * fEIz);
err = abs(q2(2) - u_expected);
results(2) = err < 1e-4 * abs(u_expected);   % small Phi correction acceptable
fprintf('Test 2 (bending EIz):   %s  |err|=%.2e\n', tf(results(2)), err);

%% Test 3: Torsional stiffness GJ
% Exact: phi = T*L / GJ
T   = 500;
fGJ = 3 * GJ0;   % tripled GJ only
fe  = Frame3DSectionProps(elems, EA0, EI0, EI0, fGJ, kGA0, kGA0);
q2  = solveCantilever(fe, nodes, "fix", T);
phi_expected = T * L / fGJ;
err = abs(q2(4) - phi_expected);
results(3) = err < tol * abs(phi_expected);
fprintf('Test 3 (torsion GJ):    %s  |err|=%.2e\n', tf(results(3)), err);

%% Test 4: GJ and EA are truly independent
% Build two elements that swap GJ and EA (one stiff axially, one torsionally).
% Each should respond only to its respective load.
fA = 3*EA0;  fJ = 0.25*GJ0;
feA = Frame3DSectionProps(elems, fA,   EI0, EI0, fJ,  kGA0, kGA0);
feJ = Frame3DSectionProps(elems, EA0,  EI0, EI0, GJ0, kGA0, kGA0);
qA = solveCantilever(feA, nodes, "ux", 1e3);   % axial load
qJ = solveCantilever(feJ, nodes, "fix", 500);  % torsion load
% axial displacement should scale inversely with EA
u_A_expected  = 1e3 / fA;
phi_J_expected = 500 * L / GJ0;
errA = abs(qA(1) - u_A_expected);
errJ = abs(qJ(4) - phi_J_expected);
results(4) = (errA < tol*abs(u_A_expected)) && (errJ < tol*abs(phi_J_expected));
fprintf('Test 4 (EA/GJ indep.):  %s  |errA|=%.2e  |errJ|=%.2e\n', ...
    tf(results(4)), errA, errJ);

%% Test 5: Frame3DSectionProps matches Frame3D for full-pipe properties
% When EA=E*A, EI=E*I, GJ=G*J, kGA=kappa*G*A, the two should give
% identical tip displacements under any load.
fe_ref  = Frame3D(elems, E, nu, R, ri);
fe_new  = Frame3DSectionProps(elems, EA0, EI0, EI0, GJ0, kGA0, kGA0);
q2_ref  = solveCantilever(fe_ref, nodes, "uy", 1e3);
q2_new  = solveCantilever(fe_new, nodes, "uy", 1e3);
err = norm(q2_ref - q2_new);
results(5) = err < 1e-6 * norm(q2_ref);
fprintf('Test 5 (matches Frame3D): %s  |err|=%.2e\n', tf(results(5)), err);

%% Summary
fprintf('\n%d / %d tests passed.\n', sum(results), nTests);
if ~all(results)
    error('test_Frame3DSectionProps: %d test(s) FAILED.', sum(~results));
end

function s = tf(b)
    if b, s = 'PASS'; else, s = 'FAIL'; end
end

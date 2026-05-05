% testReferenceModuleLoads
% Validates that ReferenceModuleSolidModel applies all six section resultants
% with < 2% integrated resultant error (after automatic scaling).
%
% Six pure unit load cases are tested:
%   N=1, Ty=1, Tz=1, My=1, Mz=1, Ms=1

clear; close all; clc;

% Ensure MorphoGenCAS_Arm classes take priority over sibling projects that
% also have a Mesh.m but lack addManipulatorHalfSegment3D.
scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));   % prepends → takes priority
clear classes;                   % flush any stale class cache

% ----- Model parameters --------------------------------------------------
arm = armModelDefaults("thin");
E = arm.E;
nu = arm.nu;
R = arm.R;
alpha_deg = arm.alpha_deg;

% Local diagnostic mesh override for load-resultant validation.
r = 0.13;
h = 0.15;
res_th = 2;

tol = 0.02;            % 2% resultant accuracy target

% ----- Build model -------------------------------------------------------
fprintf('Building ReferenceModuleSolidModel...\n');
mdl = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th);
fprintf('  Nodes : %d\n', size(mdl.mesh.nodes,1));
fprintf('  Elems : %d\n', size(mdl.mesh.elems,1));
fprintf('  xc=%.4f  yc=%.4f\n', mdl.xc, mdl.yc);
fprintf('  A=%.4e  Iy=%.4e  J=%.4e\n', mdl.A, mdl.Iy, mdl.J);

% ----- Visualise faces ---------------------------------------------------
figure('Name','Reference module - face assignment');
hold on; daspect([1 1 1]); axis on; view(3);
mdl.fe.plot(mdl.mesh.nodes);
ids_load = mdl.loaded_node_ids;
up_ids   = find(mdl.fixedFaceSelector.select(mdl.mesh.nodes));
scatter3(mdl.mesh.nodes(ids_load,1), mdl.mesh.nodes(ids_load,2), ...
         mdl.mesh.nodes(ids_load,3), 60, 'g', 'filled');
scatter3(mdl.mesh.nodes(up_ids,1), mdl.mesh.nodes(up_ids,2), ...
         mdl.mesh.nodes(up_ids,3), 60, 'r', 'filled');
legend('mesh','loaded face (green)','fixed face (red)','Location','best');
title('Reference module: loaded (green) vs fixed (red) faces');
xlabel('x'); ylabel('y'); zlabel('z');

% ----- Define unit load cases -------------------------------------------
cases = {
    'N=1',    struct('N',  1);
    'Ty=1',   struct('Ty', 1);
    'Tz=1',   struct('Tz', 1);
    'My=1',   struct('My', 1);
    'Mz=1',   struct('Mz', 1);
    'Ms=1',   struct('Ms', 1);
};

% ----- Run all cases and collect results ---------------------------------
nCases  = size(cases,1);
all_pass = true;

for k = 1:nCases
    name = cases{k,1};
    lc   = cases{k,2};
    fprintf('\n====== Load case: %s ======\n', name);
    mdl.applyLoadCase(lc);
    res = mdl.validateLoadApplication(lc, tol);
    if ~res.passed, all_pass = false; end

    % Visualise the applied nodal load
    figure('Name', ['Load: ' name]);
    hold on; daspect([1 1 1]); axis on; view(3);
    mdl.fe.plot(mdl.mesh.nodes);
    mdl.analysis.plotCurrentLoad();
    title(['Applied nodal load – ' name]);
    xlabel('x'); ylabel('y'); zlabel('z');
end

% ----- Summary -----------------------------------------------------------
fprintf('\n========================================\n');
if all_pass
    fprintf('ALL CASES PASSED (tol = %.1f%%)\n', tol*100);
else
    fprintf('SOME CASES FAILED – review output above\n');
end
fprintf('========================================\n');

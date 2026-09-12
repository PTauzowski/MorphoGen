%% test_frameBasedSolver.m
% Integration test script for ManipulatorModel3D.frameBasedSolver and
% associated methods (buildEdgeNodeSets).
%
% Run as a plain script: >> test_frameBasedSolver
%
% Each test prints PASS or FAIL with a short description.

%% --- Path setup ---------------------------------------------------------
baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(fullfile(baseDir, 'analysis')));
addpath(genpath(fullfile(baseDir, 'design')));
addpath(genpath(fullfile(baseDir, 'elements')));
addpath(genpath(fullfile(baseDir, 'examples', 'models')));
addpath(genpath(fullfile(baseDir, 'materials')));
addpath(genpath(fullfile(baseDir, 'math')));
addpath(genpath(fullfile(baseDir, 'mesh')));
addpath(genpath(fullfile(baseDir, 'examples', 'ManipulatorArm', 'HelperFunctions')));

%% --- Model parameters ---------------------------------------------------
E         = 210e9;
nu        = 0.3;
ls        = 0.15;
R         = 0.14;
r         = 0.13;
res       = 1;
res_th    = 1;
alpha_deg = 22.5;
betas_deg = [0 0 0 0 0 0 0];
Pz        = -100;
ShapeFn   = ShapeFunctionH8();
use_offset = 0;

fprintf('Building ManipulatorModel3D (minimal mesh: res=%d, res_th=%d)...\n', res, res_th);
model = ManipulatorModel3D(E, nu, ls, R, r, res, res_th, alpha_deg, betas_deg, ShapeFn, use_offset, Pz);
fprintf('Model built: %d solid nodes, %d frame nodes.\n', ...
    size(model.mesh.nodes,1), size(model.frameNodes,1));

nTests  = 5;
results = false(nTests, 1);

%% --- TEST 1: buildEdgeNodeSets returns non-empty sets -------------------
testName = 'TEST 1 - buildEdgeNodeSets: all frame node edge sets non-empty';
try
    edgeSets = model.buildEdgeNodeSets();
    nFN = size(model.frameNodes, 1);
    allNonEmpty = true;
    for k = 1:nFN
        if isempty(edgeSets{k})
            allNonEmpty = false;
            fprintf('  Frame node %d has empty edge set.\n', k);
        end
    end
    if allNonEmpty
        fprintf('PASS: %s\n', testName);
        results(1) = true;
    else
        fprintf('FAIL: %s  (some frame nodes have empty edge sets)\n', testName);
    end
catch ME
    fprintf('FAIL: %s\n  Error: %s\n', testName, ME.message);
end

%% --- TEST 2: Frame solution - end-effector displaces in load direction --
% The constructor applies [0 0 -Pz] to the frame tip node.
% With Pz = -100, the frame force is +100 in z (upward), so uz_tip > 0.
% In general the tip displacement must be non-zero and have the same sign
% as the applied frame force (-Pz).
testName = 'TEST 2 - frame solution: tip displaces in applied-load direction';
try
    nFE_frame = size(model.frame_mesh.elems, 1);
    model.frame_analysis.solveWeighted(ones(nFE_frame, 1));
    q_frame = model.frame_analysis.qnodal;   % nFrameNodes x 6
    iuz = model.frame_analysis.findDOFsIndices("uz");

    % Tip node is the last frame node
    uz_tip     = q_frame(end, iuz);
    % Frame load in z = -Pz  (constructor convention)
    frame_fz   = -Pz;
    isNonZero  = abs(uz_tip) > 1e-15;
    correctDir = (uz_tip * frame_fz) > 0;   % displacement and force same sign

    if isNonZero && correctDir
        fprintf('PASS: %s  (uz_tip = %.6e m, frame_fz = %.4g N)\n', ...
            testName, uz_tip, frame_fz);
        results(2) = true;
    elseif ~isNonZero
        fprintf('FAIL: %s  (uz_tip = 0, no displacement)\n', testName);
    else
        fprintf('FAIL: %s  (uz_tip = %.6e, frame_fz = %.4g; wrong direction)\n', ...
            testName, uz_tip, frame_fz);
    end
catch ME
    fprintf('FAIL: %s\n  Error: %s\n', testName, ME.message);
end

%% --- TEST 3: frameBasedSolver output sizes match model sizes ------------
testName = 'TEST 3 - frameBasedSolver: output sizes consistent';
try
    nSolidNodes  = size(model.mesh.nodes, 1);
    nFrameElems  = size(model.frame_mesh.elems, 1);
    nSolidElems  = model.analysis.getTotalElemsNumber();

    x = ones(nSolidElems, 1);
    [q_solid, frame_forces] = model.frameBasedSolver(x);

    sizeOK_solid = isequal(size(q_solid), [nSolidNodes, 3]);
    % frame_forces expected: 12 x nFrameElems
    sizeOK_forces = (size(frame_forces, 2) == nFrameElems) && (size(frame_forces, 1) == 12);

    if sizeOK_solid && sizeOK_forces
        fprintf('PASS: %s  (q_solid: %dx%d, frame_forces: %dx%d)\n', ...
            testName, size(q_solid,1), size(q_solid,2), ...
            size(frame_forces,1), size(frame_forces,2));
        results(3) = true;
    else
        fprintf('FAIL: %s\n', testName);
        fprintf('  q_solid size: %dx%d (expected %dx3)\n', ...
            size(q_solid,1), size(q_solid,2), nSolidNodes);
        fprintf('  frame_forces size: %dx%d (expected 12x%d)\n', ...
            size(frame_forces,1), size(frame_forces,2), nFrameElems);
    end
catch ME
    fprintf('FAIL: %s\n  Error: %s\n', testName, ME.message);
end

%% --- TEST 4: Kinematic consistency at cross-sections --------------------
testName = 'TEST 4 - kinematic consistency: mean edge displacement ~ frame displacement';
try
    % Re-fetch q_frame (frame was solved in TEST 2; frameBasedSolver re-solves it internally)
    q_frame_check = model.frame_analysis.qnodal;   % nFrameNodes x 6
    edgeSets_check = model.buildEdgeNodeSets();
    itr = model.frame_analysis.findDOFsIndices(["ux","uy","uz"]);

    nFN = size(model.frameNodes, 1);
    tol_rel = 0.50;   % 50 % relative tolerance (nodes span the cross-section)
    passCount = 0;
    skipCount = 0;

    for k = 1:nFN
        inodes = edgeSets_check{k};
        if isempty(inodes)
            skipCount = skipCount + 1;
            continue;
        end
        u_mean   = mean(model.qnodal_solid(inodes, :), 1);   % 1x3
        u_frame  = q_frame_check(k, itr);                     % 1x3

        ref = max(norm(u_frame), 1e-12);
        err = norm(u_mean - u_frame) / ref;

        if err <= tol_rel
            passCount = passCount + 1;
        else
            fprintf('  Frame node %d: err=%.3f (tol=%.2f), u_frame=%s, u_mean=%s\n', ...
                k, err, tol_rel, mat2str(u_frame,4), mat2str(u_mean,4));
        end
    end

    nChecked = nFN - skipCount;
    if nChecked > 0 && passCount == nChecked
        fprintf('PASS: %s  (%d/%d nodes checked)\n', testName, passCount, nChecked);
        results(4) = true;
    elseif nChecked == 0
        fprintf('PASS: %s  (no non-empty edge sets to check, skipped)\n', testName);
        results(4) = true;
    else
        fprintf('FAIL: %s  (%d/%d cross-sections passed tolerance)\n', ...
            testName, passCount, nChecked);
    end
catch ME
    fprintf('FAIL: %s\n  Error: %s\n', testName, ME.message);
end

%% --- TEST 5: frameBasedSolver does not mutate analysis.supports ---------
testName = 'TEST 5 - no mutation: analysis.supports unchanged after frameBasedSolver';
try
    supports_before = model.analysis.supports;

    nSolidElems = model.analysis.getTotalElemsNumber();
    model.frameBasedSolver(ones(nSolidElems, 1));

    supports_after = model.analysis.supports;

    if isequal(supports_before, supports_after)
        fprintf('PASS: %s\n', testName);
        results(5) = true;
    else
        nDiff = sum(supports_before(:) ~= supports_after(:));
        fprintf('FAIL: %s  (%d entries differ)\n', testName, nDiff);
    end
catch ME
    fprintf('FAIL: %s\n  Error: %s\n', testName, ME.message);
end

%% --- Summary ------------------------------------------------------------
fprintf('\n=== Test Summary: %d / %d passed ===\n', sum(results), nTests);
for t = 1:nTests
    if results(t)
        status = 'PASS';
    else
        status = 'FAIL';
    end
    fprintf('  Test %d: %s\n', t, status);
end

if all(results)
    fprintf('\nAll tests PASSED.\n');
else
    fprintf('\nSome tests FAILED. See output above for details.\n');
end

%% test_rotationAwareElementMapping.m
% Verifies that use_offset=false + nCircDiv produces rotation-invariant element
% numbering across arm configurations with different beta angles.
%
% The core requirement: for multi-configuration topology optimisation to be
% physically correct, design variable x(e) must represent the same local
% circumferential region of the segment in every load configuration. This holds
% iff every configuration's mesh is generated in the local segment frame (no
% phase trick), i.e. use_offset=false with resCirc snapped to nCircDiv multiples.
%
% Tests:
%   1. resCirc and halfSegmentNelems are identical across all 6 standard configs.
%   2. mesh.elems (global connectivity) is identical across configs — element
%      index e refers to the same topological position in every model.
%   3. The assertion inside ManipulatorModel3D fires when use_offset=false is
%      combined with nCircDiv=1 (missing snapping → non-conforming junctions).
%
% Run as: >> test_rotationAwareElementMapping

%% --- Path setup -----------------------------------------------------------
baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(fullfile(baseDir, 'analysis')));
addpath(genpath(fullfile(baseDir, 'design')));
addpath(genpath(fullfile(baseDir, 'elements')));
addpath(genpath(fullfile(baseDir, 'examples', 'models')));
addpath(genpath(fullfile(baseDir, 'materials')));
addpath(genpath(fullfile(baseDir, 'math')));
addpath(genpath(fullfile(baseDir, 'mesh')));
addpath(genpath(fullfile(baseDir, 'examples', 'ManipulatorArm', 'HelperFunctions')));

nPass = 0;
nFail = 0;

function reportTest(name, passed)
    if passed
        fprintf('  PASS  %s\n', name);
    else
        fprintf('  FAIL  %s\n', name);
    end
end

%% --- Shared parameters (thin preset, low resolution for speed) ------------
arm       = armModelDefaults("thin");
configs   = armLoadConfigs("six");
nConfigs  = numel(configs);

% Use a coarser mesh so building 6 models is fast.
res_th_fast = 1;
res_fast    = 5;

fprintf('Building %d models with use_offset=false, nCircDiv=%d...\n', nConfigs, arm.nCircDiv);
models = cell(nConfigs, 1);
for k = 1:nConfigs
    models{k} = ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
        res_fast, res_th_fast, arm.alpha, configs{k}.betas, arm.ShapeFn, ...
        false, arm.Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
end
fprintf('  Done.\n\n');

%% --- Test 1: resCirc is identical across all configurations ---------------
fprintf('Test 1: resCirc identical across configurations\n');
refCirc = models{1}.resCirc;
ok = true;
for k = 2:nConfigs
    if models{k}.resCirc ~= refCirc
        fprintf('  config %d (%s): resCirc=%d, expected %d\n', ...
            k, configs{k}.name, models{k}.resCirc, refCirc);
        ok = false;
    end
end
reportTest(sprintf('resCirc=%d for all %d configs', refCirc, nConfigs), ok);
if ok; nPass = nPass + 1; else; nFail = nFail + 1; end

%% --- Test 2: halfSegmentNelems identical -----------------------------------
fprintf('\nTest 2: halfSegmentNelems identical across configurations\n');
refH = models{1}.halfSegmentNelems;
ok = true;
for k = 2:nConfigs
    if models{k}.halfSegmentNelems ~= refH
        fprintf('  config %d (%s): H=%d, expected %d\n', ...
            k, configs{k}.name, models{k}.halfSegmentNelems, refH);
        ok = false;
    end
end
reportTest(sprintf('halfSegmentNelems=%d for all %d configs', refH, nConfigs), ok);
if ok; nPass = nPass + 1; else; nFail = nFail + 1; end

%% --- Test 3: mesh.elems (connectivity) identical --------------------------
fprintf('\nTest 3: mesh.elems connectivity identical across configurations\n');
refElems = models{1}.mesh.elems;
ok = true;
for k = 2:nConfigs
    if ~isequal(models{k}.mesh.elems, refElems)
        nDiff = nnz(models{k}.mesh.elems ~= refElems);
        fprintf('  config %d (%s): elems differ in %d entries\n', ...
            k, configs{k}.name, nDiff);
        ok = false;
    end
end
reportTest(sprintf('mesh.elems identical for all %d configs', nConfigs), ok);
if ok; nPass = nPass + 1; else; nFail = nFail + 1; end

%% --- Test 4: element count matches across configurations ------------------
fprintf('\nTest 4: total element count identical across configurations\n');
refNElems = models{1}.analysis.getTotalElemsNumber();
ok = true;
for k = 2:nConfigs
    nk = models{k}.analysis.getTotalElemsNumber();
    if nk ~= refNElems
        fprintf('  config %d (%s): nElems=%d, expected %d\n', ...
            k, configs{k}.name, nk, refNElems);
        ok = false;
    end
end
reportTest(sprintf('nElems=%d for all %d configs', refNElems, nConfigs), ok);
if ok; nPass = nPass + 1; else; nFail = nFail + 1; end

%% --- Test 5: interface node sharing matches assertLinkedArmLayoutCompatible
fprintf('\nTest 5: adjacent half-segment interfaces are conforming\n');
ok = true;
try
    for k = 1:nConfigs
        assertLinkedArmLayoutCompatible(models{k}, models{1}, configs{k}.name);
    end
catch ME
    fprintf('  Error: %s\n', ME.message);
    ok = false;
end
reportTest('All configs pass assertLinkedArmLayoutCompatible', ok);
if ok; nPass = nPass + 1; else; nFail = nFail + 1; end

%% --- Test 6: use_offset=false + nCircDiv=1 fires the guard ----------------
fprintf('\nTest 6: ManipulatorModel3D asserts when use_offset=false, nCircDiv=1\n');
ok = false;
try
    ManipulatorModel3D(arm.E, arm.nu, arm.h_seg, arm.R, arm.r, ...
        res_fast, res_th_fast, arm.alpha, configs{1}.betas, arm.ShapeFn, ...
        false, arm.Pz, arm.constEndRing, arm.constMiddleRing, 1);
    % Should not reach here
catch ME
    ok = contains(ME.message, 'nCircDiv');
    if ~ok
        fprintf('  Unexpected error: %s\n', ME.message);
    end
end
reportTest('Guard fires with nCircDiv=1 and use_offset=false', ok);
if ok; nPass = nPass + 1; else; nFail = nFail + 1; end

%% --- Summary ---------------------------------------------------------------
fprintf('\n========================================\n');
fprintf('Results: %d passed, %d failed\n', nPass, nFail);
if nFail == 0
    fprintf('ALL TESTS PASSED\n');
else
    fprintf('SOME TESTS FAILED\n');
end
fprintf('========================================\n');

if nFail > 0
    error('test_rotationAwareElementMapping: %d test(s) failed.', nFail);
end

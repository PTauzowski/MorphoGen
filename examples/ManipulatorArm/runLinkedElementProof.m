% runLinkedElementProof
% Runner script for plotLinkedElement.
%
% Builds one full-arm ManipulatorModel3D (single configuration, the
% straight-bending case is cheapest to construct) and then generates
% three figures that prove the segmentToArm element-linkage:
%
%   Figure 1 — randomly chosen element (visual surprise / sanity check)
%   Figure 2 — element near the circumferential mid-point of the tube
%   Figure 3 — element from the inner-radial layer (next to void interior)
%
% All figures are saved to results/linkedElementProof/.
%
% The key assertion being demonstrated:
%   Changing rho(e_ref) changes EXACTLY nCopies elements in x_arm,
%   which are the nArms red (2a) and nArms blue (2b) highlighted bands.
%   Every configuration with a 45° joint (max_torsion, min_torsion) will
%   produce the same connectivity — proving that joint rotations do NOT
%   break element correspondence across segments.

clear; close all; clc;
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

resultRoot = fullfile(scriptDir, 'results', 'linkedElementProof');
if ~exist(resultRoot, 'dir')
    mkdir(resultRoot);
end

rng(42, 'twister');

%% ---- Arm geometry ----------------------------------------------------------
arm    = armModelDefaults("thin");
E      = arm.E;
nu     = arm.nu;
R      = arm.R;
r      = arm.r;
h_seg  = arm.h_seg;
alpha  = arm.alpha;
res    = arm.res;
res_th = arm.res_th;
Pz     = arm.Pz;
ShapeFn = arm.ShapeFn;

%% ---- Build one model per configuration ------------------------------------
% We use the "six" set: max_torsion includes a 45-degree joint, which is
% the hardest case for element-alignment across segments.
configs = armLoadConfigs("six");

fprintf('Building %d arm configurations...\n', numel(configs));
models = cell(numel(configs), 1);
for k = 1:numel(configs)
    cfg = configs{k};
    fprintf('  [%d/%d] %s  betas=%s\n', k, numel(configs), cfg.label, mat2str(cfg.betas));
    models{k} = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
        cfg.betas, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);
end

H      = models{1}.halfSegmentNelems;
nElems = models{1}.analysis.getTotalElemsNumber();
nCopies = nElems / H;
nArms   = nCopies / 2;
resCirc = models{1}.resCirc;

fprintf('\nArm statistics:\n');
fprintf('  H (half-segment elements) : %d\n', H);
fprintf('  nElems (full arm)         : %d\n', nElems);
fprintf('  nCopies                   : %d\n', nCopies);
fprintf('  nArms                     : %d\n', nArms);
fprintf('  resCirc (snapped)         : %d  (nCircDiv=%d)\n', resCirc, arm.nCircDiv);

%% ---- Verify integer circumferential shifts for all configurations ----------
fprintf('\nCircumferential shift check (must be integer for exact local-frame linking):\n');
cfg1 = configs{1};
betas1_rad = cfg1.betas * pi / 180;
cumPhi = cumsum(betas1_rad);
shifts = cumPhi * resCirc / (2*pi);
allInteger = all(abs(shifts - round(shifts)) < 1e-9);
for s = 1:numel(shifts)
    fprintf('  segment %d: cumPhi=%.1f deg  shift=%.4f  integer=%d\n', ...
        s, cumPhi(s)*180/pi, shifts(s), abs(shifts(s)-round(shifts(s)))<1e-9);
end
assert(allInteger, 'Non-integer circumferential shifts detected — increase nCircDiv or adjust joint angles.');

%% ---- Verify that ALL configurations share the same linked layout ----------
fprintf('\nLinked-layout check across all %d configurations:\n', numel(configs));
refH       = models{1}.halfSegmentNelems;
refNElems  = models{1}.analysis.getTotalElemsNumber();
refNCopies = refNElems / refH;
refResCirc = models{1}.resCirc;
refResLen  = models{1}.resLen;
refResTh   = models{1}.resTh;

for k = 1:numel(configs)
    m = models{k};
    Hk = m.halfSegmentNelems;
    nElemsK = m.analysis.getTotalElemsNumber();
    nCopiesK = nElemsK / Hk;

    sameLayout = Hk == refH && nElemsK == refNElems && ...
        nCopiesK == refNCopies && ...
        m.resCirc == refResCirc && m.resLen == refResLen && m.resTh == refResTh;

    fprintf('  %s  (betas=%s)  same linked layout: %d\n', ...
        configs{k}.label, mat2str(configs{k}.betas), sameLayout);
    assert(sameLayout, 'Linked layout mismatch in configuration %s!', configs{k}.name);

    expectedShared = 2 * m.resCirc * m.resTh;
    for c = 1:nCopiesK-1
        e1 = (c-1)*Hk + (1:Hk);
        e2 = c*Hk + (1:Hk);
        shared = intersect(unique(m.mesh.elems(e1,:)), unique(m.mesh.elems(e2,:)));
        assert(numel(shared) == expectedShared, ...
            'Interface %d-%d in %s has %d shared nodes, expected %d.', ...
            c, c+1, configs{k}.name, numel(shared), expectedShared);
    end
end
fprintf('  PASS: all configurations share linked layout and conforming interfaces.\n');

%% ---- Choose three illustrative reference elements -------------------------
% 1. Random — for the "surprise" demo
e_random = randi(H);

% 2. Circumferential mid-point of the outer radial layer, axial mid
%    Elements are ordered (ix, iy, iz): ix=radial, iy=circumferential, iz=axial
%    Uses the actual snapped mesh resolution stored on the model.
%    ix=1 only, iy=resCirc/2, iz=resLen/2
modelForIndexing = models{1};
resTh   = modelForIndexing.resTh;
resCirc = modelForIndexing.resCirc;
resLen  = modelForIndexing.resLen;
e_circ  = sub2elemIndex(resTh, resCirc, resLen, 1, round(resCirc/2), round(resLen/2));

% 3. Same circumferential position but axial-end (closest to the junction ring)
e_end   = sub2elemIndex(resTh, resCirc, resLen, 1, round(resCirc/2), 1);

elemChoices = {e_random, 'random'; ...
               e_circ,   'circ_mid'; ...
               e_end,    'end_ring'};

fprintf('\nSelected reference elements:\n');
fprintf('  e_random = %d\n', e_random);
fprintf('  e_circ   = %d  (outer layer, circumferential mid, axial mid)\n', e_circ);
fprintf('  e_end    = %d  (outer layer, circumferential mid, axial end)\n', e_end);

%% ---- Generate and save figures --------------------------------------------
% Use the first model (max_bending — straight, easy to read visually)
% but any model would look identical because connectivity is the same.
model = models{1};
baseCfg = configs{1};

fprintf('\nGenerating linked-element figures (model: %s, betas=%s)...\n', ...
    baseCfg.label, mat2str(baseCfg.betas));

for i = 1:size(elemChoices, 1)
    e_ref  = elemChoices{i, 1};
    label  = elemChoices{i, 2};

    fprintf('\n--- Element %d (%s) ---\n', e_ref, label);
    fig = plotLinkedElement(model, e_ref, resultRoot, baseCfg.label, baseCfg.betas);
    set(fig, 'Name', sprintf('%s linked element %d (%s)', baseCfg.name, e_ref, label));
    % rename saved files to include the label
    oldPng = fullfile(resultRoot, sprintf('linked_element_%d.png', e_ref));
    oldFig = fullfile(resultRoot, sprintf('linked_element_%d.fig', e_ref));
    newPng = fullfile(resultRoot, sprintf('linked_element_%s_%s_%d.png', ...
        baseCfg.name, label, e_ref));
    newFig = fullfile(resultRoot, sprintf('linked_element_%s_%s_%d.fig', ...
        baseCfg.name, label, e_ref));
    if exist(oldPng, 'file'), movefile(oldPng, newPng); end
    if exist(oldFig, 'file'), movefile(oldFig, newFig); end
end

%% ---- Connectivity stress-test: generate figure for max_torsion (45 deg) --
% This is the configuration the colleague questioned — 45-degree joints.
% The figure must look identical structurally (same element bands).
fprintf('\n--- Connectivity stress-test: max_torsion (45-degree joints) ---\n');
idx_torsion = find(cellfun(@(c) strcmp(c.name, 'max_torsion'), configs), 1);
if ~isempty(idx_torsion)
    torsionCfg = configs{idx_torsion};
    fig_torsion = plotLinkedElement(models{idx_torsion}, e_circ, resultRoot, ...
        torsionCfg.label, torsionCfg.betas);
    set(fig_torsion, 'Name', sprintf('%s e_ref=%d', torsionCfg.name, e_circ));
    oldPng = fullfile(resultRoot, sprintf('linked_element_%d.png', e_circ));
    oldFig = fullfile(resultRoot, sprintf('linked_element_%d.fig', e_circ));
    if exist(oldPng, 'file')
        movefile(oldPng, fullfile(resultRoot, sprintf('linked_element_%s_%d.png', ...
            torsionCfg.name, e_circ)));
    end
    if exist(oldFig, 'file')
        movefile(oldFig, fullfile(resultRoot, sprintf('linked_element_%s_%d.fig', ...
            torsionCfg.name, e_circ)));
    end
    fprintf('  Saved torsion figure (betas=%s).\n', mat2str(torsionCfg.betas));
end

fprintf('\nAll figures saved to %s\n', resultRoot);
fprintf('Each red band = one 2a (normal) copy of rho(%d).\n', e_circ);
fprintf('Each blue band = one 2b (flipped) copy at mirror element H+1-%d=%d.\n', ...
    e_circ, H + 1 - e_circ);

%% ---- Local helper ----------------------------------------------------------
function idx = sub2elemIndex(resTh, resCirc, resLen, ix, iy, iz)
% Convert (ix, iy, iz) subscript (1-based, ix fastest) to 1-based element index.
%   ix in [1..resTh], iy in [1..resCirc], iz in [1..resLen]
    ix = min(max(ix, 1), resTh);
    iy = min(max(iy, 1), resCirc);
    iz = min(max(iz, 1), resLen);
    idx = (iz - 1) * resTh * resCirc + (iy - 1) * resTh + ix;
end

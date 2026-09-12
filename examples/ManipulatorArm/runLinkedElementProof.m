% runLinkedElementProof
% Runner script for plotLinkedElement.
%
% Builds one full-arm ManipulatorModel3D per configuration from
% armLoadConfigs("sixPlusTension") and generates two families of figures:
%
%   Section 1 — Rotation-aware mesh visualization (NEW)
%     One figure per configuration showing the FULL OPAQUE mesh (gray faces,
%     black edges) with the "beta=0 stripe" highlighted in red.  The stripe
%     is the set of all elements at circumferential local index iy=1 (angle=0
%     in each segment's own local frame) across every half-segment copy.
%     Because use_offset=false generates every segment in its local frame,
%     the stripe runs parallel to each segment axis but steps at each joint
%     by the cumulative beta rotation — directly visualising segment rotation
%     emulation.  Saved as: results/linkedElementProof/rotation_proof_NN_<name>.{fig,png}
%
%   Section 2 — Linked element proof (existing)
%     Three figures that prove the segmentToArm element-linkage by
%     colouring all full-arm copies of one reference half-segment element.
%     The key assertion: changing rho(e_ref) changes EXACTLY nCopies elements
%     in x_arm, which are the nArms red (2a) and nArms blue (2b) highlighted
%     bands.  Every configuration with a 45-degree joint (max_torsion) will
%     produce the same linked layout — proving joint rotations do NOT break
%     element correspondence across segments.

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
configs = armLoadConfigs("sixPlusTension");

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
resTh_m = models{1}.resTh;
resLen_m = models{1}.resLen;

fprintf('\nArm statistics:\n');
fprintf('  H (half-segment elements) : %d\n', H);
fprintf('  nElems (full arm)         : %d\n', nElems);
fprintf('  nCopies                   : %d\n', nCopies);
fprintf('  nArms                     : %d\n', nArms);
fprintf('  resCirc (snapped)         : %d  (nCircDiv=%d)\n', resCirc, arm.nCircDiv);
fprintf('  resTh                     : %d\n', resTh_m);
fprintf('  resLen                    : %d\n', resLen_m);

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
    assertLinkedArmLayoutCompatible(m, models{1}, configs{k}.name);
end
fprintf('  PASS: all configurations share linked layout and conforming interfaces.\n');

%% ========================================================================
%  Section 1 — Rotation-aware mesh visualization
%  One figure per configuration: full opaque mesh (gray) + red beta=0 stripe.
%
%  "Beta=0 stripe": elements at circumferential local index iy=1 (angle=0 in
%  each segment's local frame) across ALL half-segment copies.  This forms a
%  band running parallel to each segment axis.  For rotated joints the band
%  steps by the joint beta — directly showing the rotation emulation.
% =========================================================================
fprintf('\n=== Section 1: Rotation-aware mesh visualisation ===\n');

% Local indices within one half-segment for iy=1 (angle=0 in local frame).
% Element ordering: idx = (iz-1)*resTh*resCirc + (iy-1)*resTh + ix
%   iy=1 → idx = (iz-1)*resTh*resCirc + ix,  ix=1..resTh, iz=1..resLen
localRedIds = reshape( ...
    bsxfun(@plus, (0:resLen_m-1)' * resTh_m * resCirc, 1:resTh_m), ...
    [], 1);   % [resTh*resLen  x  1]
fprintf('  Red stripe: %d elements per half-segment copy (iy=1, all radial, all axial)\n', ...
    numel(localRedIds));

for k = 1:numel(configs)
    cfg   = configs{k};
    model = models{k};
    nE    = size(model.mesh.elems, 1);
    nCop  = nE / H;

    % Global element ids for the beta=0 stripe across every half-segment copy.
    % For copy c (0-based): global id = c*H + localRedId
    globalRedIds = reshape( ...
        bsxfun(@plus, (0:nCop-1) * H, localRedIds), ...
        [], 1);   % [resTh*resLen*nCop  x  1]
    maskRed = false(nE, 1);
    maskRed(globalRedIds) = true;

    fig = figure('Color', 'white', 'Units', 'normalized', ...
        'Position', [0.05 0.05 0.55 0.75]);
    hold on; axis off; daspect([1 1 1]);
    view([30 20]);
    light('Position', [-1 -2 5], 'Style', 'infinite');
    light('Position', [ 1  1 3], 'Style', 'local');
    lighting flat;
    material dull;

    % Single patch: exterior faces only, per-face flat color (gray / red stripe).
    % One patch object = fast interactive rotation.
    plotMeshTwoColor(model.mesh.nodes, model.mesh.elems, model.fe.sf.fcontours, ...
        maskRed, [0.65 0.65 0.65], [0.85 0.15 0.10]);

    betaStr = mat2str(cfg.betas);
    title({cfg.label, sprintf('betas = %s', betaStr)}, ...
        'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');

    % Save
    stem = fullfile(resultRoot, sprintf('rotation_proof_%02d_%s', k, cfg.name));
    exportgraphics(fig, [stem '.png'], 'Resolution', 200);
    savefig(fig, [stem '.fig']);
    fprintf('  [%d/%d] %s  -->  %s\n', k, numel(configs), cfg.label, stem);
    close(fig);
end
fprintf('Section 1 complete.\n');

%% ========================================================================
%  Section 2 — Linked element proof (plotLinkedElement figures)
% =========================================================================
fprintf('\n=== Section 2: Linked element proof ===\n');

%% ---- Choose three illustrative reference elements -------------------------
% 1. Random — for the "surprise" demo
e_random = randi(H);

% 2. Circumferential mid-point of the outer radial layer, axial mid
modelForIndexing = models{1};
e_circ  = sub2elemIndex(resTh_m, resCirc, resLen_m, 1, round(resCirc/2), round(resLen_m/2));

% 3. Same circumferential position but axial-end (closest to the junction ring)
e_end   = sub2elemIndex(resTh_m, resCirc, resLen_m, 1, round(resCirc/2), 1);

elemChoices = {e_random, 'random'; ...
               e_circ,   'circ_mid'; ...
               e_end,    'end_ring'};

fprintf('\nSelected reference elements:\n');
fprintf('  e_random = %d\n', e_random);
fprintf('  e_circ   = %d  (outer layer, circumferential mid, axial mid)\n', e_circ);
fprintf('  e_end    = %d  (outer layer, circumferential mid, axial end)\n', e_end);

%% ---- Generate and save figures --------------------------------------------
model   = models{1};
baseCfg = configs{1};

fprintf('\nGenerating linked-element figures (model: %s, betas=%s)...\n', ...
    baseCfg.label, mat2str(baseCfg.betas));

for i = 1:size(elemChoices, 1)
    e_ref  = elemChoices{i, 1};
    label  = elemChoices{i, 2};

    fprintf('\n--- Element %d (%s) ---\n', e_ref, label);
    fig = plotLinkedElement(model, e_ref, resultRoot, baseCfg.label, baseCfg.betas);
    set(fig, 'Name', sprintf('%s linked element %d (%s)', baseCfg.name, e_ref, label));
    oldPng = fullfile(resultRoot, sprintf('linked_element_%d.png', e_ref));
    oldFig = fullfile(resultRoot, sprintf('linked_element_%d.fig', e_ref));
    newPng = fullfile(resultRoot, sprintf('linked_element_%s_%s_%d.png', ...
        baseCfg.name, label, e_ref));
    newFig = fullfile(resultRoot, sprintf('linked_element_%s_%s_%d.fig', ...
        baseCfg.name, label, e_ref));
    if exist(oldPng, 'file'), movefile(oldPng, newPng); end
    if exist(oldFig, 'file'), movefile(oldFig, newFig); end
end

%% ---- Connectivity stress-test: max_torsion (45 deg) -----------------------
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

fprintf('\nSection 2 complete.\n');
fprintf('\nAll outputs saved to %s\n', resultRoot);
fprintf('Section 1: rotation_proof_NN_<name>.{fig,png}  — one per configuration\n');
fprintf('Section 2: linked_element_*.{fig,png}           — linked element proof\n');
fprintf('\nRed stripe (Section 1): iy=1 across all %d half-segment copies.\n', nCopies);
fprintf('Each copy contributes %d red elements (%d radial × %d axial).\n', ...
    numel(localRedIds), resTh_m, resLen_m);
fprintf('Each red band = one 2a (normal) copy of rho(%d).\n', e_circ);
fprintf('Each blue band = one 2b (flipped) copy at mirror element H+1-%d=%d.\n', ...
    e_circ, H + 1 - e_circ);

%% ---- Local helpers ---------------------------------------------------------
function idx = sub2elemIndex(resTh, resCirc, resLen, ix, iy, iz)
% Convert (ix, iy, iz) subscript (1-based, ix fastest) to 1-based element index.
%   ix in [1..resTh], iy in [1..resCirc], iz in [1..resLen]
    ix = min(max(ix, 1), resTh);
    iy = min(max(iy, 1), resCirc);
    iz = min(max(iz, 1), resLen);
    idx = (iz - 1) * resTh * resCirc + (iy - 1) * resTh + ix;
end

function plotMeshTwoColor(nodes, elems, fcon, maskB, colorA, colorB)
% Render mesh exterior as a single patch for fast interactive rotation.
% Elements with maskB=true get colorB; all others get colorA.
% fcon is sf.fcontours: [nVertsPerFace x nFacesPerElem] (columns = faces).
% Face list built with the same reshape formula as plotSolidSelected.
    nE   = size(elems, 1);
    nVpf = size(fcon, 1);   % rows  = vertices per face (4 for hex8)
    nFpe = size(fcon, 2);   % cols  = faces per element (6 for hex8)

    % allFaces: [nE*nFpe x nVpf] — faces 1..nFpe of elem 1, then elem 2, …
    allFaces  = reshape(elems(:, fcon)', nVpf, nFpe * nE)';
    ownerElem = repelem((1:nE)', nFpe);   % element that owns each row

    % Keep only exterior faces — those that appear exactly once.
    [~, ~, ic] = unique(sort(allFaces, 2), 'rows', 'stable');
    cnt  = accumarray(ic, 1);
    isExt    = cnt(ic) == 1;
    extFaces = allFaces(isExt, :);
    extOwner = ownerElem(isExt);

    % Per-face flat color.
    fc = repmat(colorA, size(extFaces, 1), 1);
    fc(maskB(extOwner), :) = repmat(colorB, nnz(maskB(extOwner)), 1);

    patch('Vertices', nodes, 'Faces', extFaces, ...
        'FaceVertexCData', fc, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 1.0);
end

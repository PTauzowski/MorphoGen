% stage3A_fullArmVerification.m
% Stage 3A: Full-arm verification of reference-module SIMP designs.
%
% Loads two optimised reference-module topologies (bending-only and
% multi-load shell), expands each to the full Arm-Z arm via segmentToArm,
% then runs frameBasedSolver for five representative configurations.
%
% Designs:
%   A  My_only (vf=0.40, Rmin=0.5t)
%   B  My+Mz+Ms+Ty+Tz (vf=0.40, Rmin=0.5t, p=4)  -- multi-load shell
%   C  Full solid (rho=1 everywhere, reference baseline)
%
% Configurations (7-joint arm, betas in degrees):
%   straight | max_bending | max_torsion | max_shear | multi_rep
%
% Pass criteria:
%   - solver converges for all design/config combinations
%   - shell design (B) <= bending design (A) in maxHM under torsion/multi configs
%   - bending design (A) may be competitive in bending configs

clear; close all; clc;
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(0, 'twister');

%% ---- Model geometry (must match reference-module Stage 2) ---------------
E      = 2.0e9;     % Pa  (PLA, same as Stage 2)
nu     = 0.35;
R      = 0.14;      % m outer radius
r      = 0.08;      % m inner radius
h_seg  = 0.25;      % m segment length
alpha  = 22.5;      % deg inclination
res_th = 4;         % radial layers (matches reference module: 4012 elems/half-seg)
Pz     = 100;       % N  tip force applied to frame
ShapeFn = ShapeFunctionL8();

%% ---- Load module designs -------------------------------------------------
sweepRoot = fullfile(scriptDir, 'results', 'referenceModuleSIMP_sweep');

fprintf('Loading reference-module designs from Stage 2C sweep...\n');

% Design A: bending-only (My only, vf=0.40, Rmin=0.5t)
dA = load(fullfile(sweepRoot, 'My_only_vf040_rmin050', 'result.mat'), 'rho_opt');
rhoA = dA.rho_opt;

% Design B: multi-load shell (My+Mz+Ms+Ty+Tz, vf=0.40, Rmin=0.5t, p=4)
dB = load(fullfile(sweepRoot, ...
    'My_Mz_Ms_Ty_Tz_vf040_rmin050_p04_vf040_rmin050_p04', 'result.mat'), 'rho_opt');
rhoB = dB.rho_opt;

nHalfSegElems = numel(rhoA);
assert(numel(rhoB) == nHalfSegElems, ...
    'Designs A and B have different half-segment element counts (%d vs %d).', ...
    numel(rhoA), numel(rhoB));
fprintf('  Half-segment element count: %d\n', nHalfSegElems);

designs = {
    struct('name', 'A_My_only',         'label', 'My only',        'rhoSeg', rhoA, 'color', [0.20 0.55 1.00]);
    struct('name', 'B_multiload_shell',  'label', 'Multi-load shell','rhoSeg', rhoB, 'color', [1.00 0.40 0.15]);
    struct('name', 'C_full_solid',       'label', 'Full solid',     'rhoSeg', [],   'color', [0.35 0.75 0.35]);
};
nDes = numel(designs);

%% ---- Arm configurations (betas in degrees, 7 joints) --------------------
% Taken from coupledExamples.m named configurations.
configs = {
    struct('name', 'straight',    'label', 'Straight',      'betas', [0   0   0   0   0   0   0]);
    struct('name', 'max_bending', 'label', 'Max M_z',       'betas', [0   0   0 180 180 180 180]);
    struct('name', 'max_torsion', 'label', 'Max M_s',       'betas', [0  45  45  45 270 180 180]);
    struct('name', 'max_shear',   'label', 'Max T_y',       'betas', [0   0 180   0 180 180 180]);
    struct('name', 'multi_rep',   'label', 'Min N / multi', 'betas', [0 180 180 180 180 180 180]);
};
nConf = numel(configs);

%% ---- Result directory ---------------------------------------------------
resultRoot = fullfile(scriptDir, 'results', 'stage3A_fullArmVerification');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

fprintf('\nFull-arm model will have ~%d elements per configuration.\n', ...
    nHalfSegElems * 12);   % 7-joint arm → 12 half-segments

%% ---- Main sweep ---------------------------------------------------------
% Preallocate metric arrays: [nConf × nDes]
tipDisp_mm  = nan(nConf, nDes);
maxDisp_mm  = nan(nConf, nDes);
maxHM_MPa   = nan(nConf, nDes);
volFracs    = nan(nConf, nDes);
statuses    = strings(nConf, nDes);

for ci = 1:nConf
    cfg = configs{ci};
    fprintf('\n=== Config %d/%d: %s  betas=%s ===\n', ...
        ci, nConf, cfg.label, mat2str(cfg.betas));

    for di = 1:nDes
        des = designs{di};
        fprintf('  [%s] ... ', des.label);
        tic;
        try
            %% Build full-arm model for this configuration
            model = ManipulatorModel3D(E, nu, h_seg, R, r, 15, res_th, alpha, ...
                cfg.betas, ShapeFn, true, Pz);
            nArmElems = model.analysis.getTotalElemsNumber();

            %% Expand module design to full arm
            if isempty(des.rhoSeg)
                x_arm = ones(nArmElems, 1);           % full solid baseline
            else
                x_arm = model.segmentToArm(des.rhoSeg);
                assert(numel(x_arm) == nArmElems, ...
                    'segmentToArm returned %d elements, arm model has %d.', ...
                    numel(x_arm), nArmElems);
            end

            %% Solve: frame → kinematic BCs → static condensation
            % frameBasedSolver also calls computeElementResults internally.
            [~, ~] = model.frameBasedSolver(x_arm);

            %% Compute metrics
            % Displacements [nNodes × 3] in metres
            qn = model.analysis.qnodal;
            dispNorm_m = sqrt(sum(qn.^2, 2));
            maxDisp_mm(ci, di)  = max(dispNorm_m) * 1e3;

            tipIds = find(model.loadSurfaceNodes);
            tipDisp_mm(ci, di)  = norm(mean(qn(tipIds, :), 1)) * 1e3;

            % Huber-Mises stress: gp.all is [nResults × nElems × nGP]
            % index 13 = sHM, restrict to active elements (rho > 0.5)
            activeElems = x_arm > 0.5;
            if ~any(activeElems)
                activeElems = true(nArmElems, 1);
            end
            hmGP = model.fe.results.gp.all(13, activeElems, :);
            maxHM_MPa(ci, di) = max(hmGP(:)) / 1e6;

            volFracs(ci, di) = mean(x_arm);

            statuses(ci, di) = "ok";
            elapsed = toc;
            fprintf('tipDisp=%.3f mm  maxHM=%.4f MPa  vf=%.3f  (%.1f s)\n', ...
                tipDisp_mm(ci, di), maxHM_MPa(ci, di), volFracs(ci, di), elapsed);

            %% Save HM stress figure for this (config, design) pair
            figPath = fullfile(resultRoot, sprintf('%s__%s', cfg.name, des.name));
            saveStressFigure(model, x_arm, activeElems, cfg.label, des.label, ...
                maxHM_MPa(ci, di), tipDisp_mm(ci, di), figPath);

        catch ME
            statuses(ci, di) = "failed";
            fprintf('FAILED: %s\n', ME.message);
            save(fullfile(resultRoot, ...
                sprintf('failed_%s_%s.mat', cfg.name, des.name)), 'ME', 'cfg', 'des');
        end
    end
end

%% ---- Summary table -------------------------------------------------------
rows = {};
for ci = 1:nConf
    for di = 1:nDes
        row.configName  = string(configs{ci}.name);
        row.configLabel = string(configs{ci}.label);
        row.designName  = string(designs{di}.name);
        row.designLabel = string(designs{di}.label);
        row.status      = statuses(ci, di);
        row.volFrac     = volFracs(ci, di);
        row.tipDisp_mm  = tipDisp_mm(ci, di);
        row.maxDisp_mm  = maxDisp_mm(ci, di);
        row.maxHM_MPa   = maxHM_MPa(ci, di);
        rows{end+1} = row; %#ok<SAGROW>
    end
end
T = struct2table(vertcat(rows{:}));
summaryFile = fullfile(resultRoot, 'summary.csv');
writetable(T, summaryFile);
fprintf('\nSummary saved to %s\n', summaryFile);
disp(T(:, {'configLabel','designLabel','volFrac','tipDisp_mm','maxHM_MPa','status'}));

%% ---- Pass/fail check ----------------------------------------------------
fprintf('\n--- Pass criteria ---\n');

% 1. All solver runs converged
allOk = all(statuses == "ok", 'all');
fprintf('  All configs/designs solved OK : %d\n', allOk);

% 2. Shell design (B) <= bending-only (A) in maxHM for torsion/multi configs
configNames = cellfun(@(c) c.name, configs, 'UniformOutput', false);
torsionIdx  = find(strcmp(configNames, 'max_torsion'));
multiIdx    = find(strcmp(configNames, 'multi_rep'));
checkCfgs   = [torsionIdx, multiIdx];
iA = 1; iB = 2;
shellBetter = all(maxHM_MPa(checkCfgs, iB) <= maxHM_MPa(checkCfgs, iA));
fprintf('  Shell design <= bending design HM (torsion+multi) : %d\n', shellBetter);

% 3. Bending-only design may be competitive in bending
bendIdx = find(strcmp(configNames, 'max_bending'));
bendCompetitive = maxHM_MPa(bendIdx, iA) <= maxHM_MPa(bendIdx, iB) * 1.5;
fprintf('  Bending design competitive in bending (within 1.5x) : %d\n', bendCompetitive);

%% ---- Comparison figures -------------------------------------------------
configLabels = cellfun(@(c) c.label, configs, 'UniformOutput', false);
designLabels = cellfun(@(d) d.label, designs, 'UniformOutput', false);
colors       = cellfun(@(d) d.color, designs, 'UniformOutput', false);

% Figure 1: tip displacement comparison
fig1 = figure('Visible', 'off', 'Name', 'tipDisp comparison', ...
    'Position', [100 100 900 420]);
bar3A_comparison(tipDisp_mm, configLabels, designLabels, colors, ...
    'Tip displacement (mm)', 'Stage 3A: Tip displacement by design and configuration');
saveas(fig1, fullfile(resultRoot, 'comparison_tipDisp.png'));
close(fig1);

% Figure 2: max HM stress comparison
fig2 = figure('Visible', 'off', 'Name', 'maxHM comparison', ...
    'Position', [100 100 900 420]);
bar3A_comparison(maxHM_MPa, configLabels, designLabels, colors, ...
    'Max HM stress (MPa)', 'Stage 3A: Max Huber-Mises stress by design and configuration');
saveas(fig2, fullfile(resultRoot, 'comparison_maxHM.png'));
close(fig2);

% Figure 3: stress ratio A/B (>1 means bending design is worse)
fig3 = figure('Visible', 'off', 'Name', 'stress ratio', 'Position', [100 100 900 420]);
ratioAB = maxHM_MPa(:, 1) ./ maxHM_MPa(:, 2);
bar(ratioAB, 'FaceColor', [0.6 0.6 0.9]);
set(gca, 'XTick', 1:nConf, 'XTickLabel', configLabels, 'XTickLabelRotation', 20);
yline(1.0, '--k', 'equal stress', 'LabelHorizontalAlignment', 'left');
ylabel('max HM ratio  A / B');
title('Stage 3A: Bending-only vs. shell design stress ratio (>1 = shell is better)');
grid on;
saveas(fig3, fullfile(resultRoot, 'stress_ratio_A_over_B.png'));
close(fig3);

fprintf('\nComparison figures saved to %s\n', resultRoot);

% =========================================================================
function saveStressFigure(model, x_arm, activeElems, cfgLabel, desLabel, ...
        maxHM_MPa, tipDisp_mm, figPath)
    fig = figure('Visible', 'off', 'Name', [cfgLabel ' - ' desLabel]);
    hold on; axis on; daspect([1 1 1]); view(45, 25);

    % Node-averaged HM stress (MPa) from element results
    hmNodal = model.fe.results.nodal.all(:, 13) / 1e6;

    elems = model.mesh.elems;
    nodes = model.mesh.nodes;
    faceP = model.fe.sf.fcontours';      % [nFacesPerElem × nodesPerFace]
    nFPE  = size(faceP, 1);
    nFPN  = size(faceP, 2);
    actIds = find(activeElems);
    nAct  = numel(actIds);

    % Vectorized face assembly: iterate only over nFPE face patterns (6 for hex8)
    nFacesTotal = nAct * nFPE;
    allFaces  = zeros(nFacesTotal, nFPN);
    allColors = zeros(nFacesTotal, 1);
    for fi = 1:nFPE
        r0 = (fi-1)*nAct + 1;  r1 = fi*nAct;
        localNids = faceP(fi, :);              % 1 × nFPN node slots in hex
        allFaces(r0:r1, :) = elems(actIds, localNids);
        faceNodeHM = hmNodal(elems(actIds, localNids));   % nAct × nFPN
        allColors(r0:r1) = mean(faceNodeHM, 2);
    end

    patch('Vertices', nodes, 'Faces', allFaces, ...
        'FaceVertexCData', allColors, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', 'FaceAlpha', 1.0);
    colormap(jet); cb = colorbar;
    cb.Label.String = 'HM stress (MPa)';
    maxC = max(allColors);
    if maxC > 0, caxis([0 maxC]); end

    % Overlay deformed arm skeleton (amplified ×50)
    scale = 50;
    dof_tr = model.frame_analysis.findDOFsIndices(["ux","uy","uz"]);
    qf     = model.frame_analysis.qnodal;    % [nFrameNodes × 6]
    fnDef  = model.frameNodes + scale * qf(:, dof_tr);
    plot3(fnDef(:,1), fnDef(:,2), fnDef(:,3), 'k-o', ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k');

    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('%s | %s\nmaxHM=%.4f MPa, tipDisp=%.3f mm (deform \times%d)', ...
        cfgLabel, desLabel, maxHM_MPa, tipDisp_mm, scale), 'FontSize', 9);
    saveas(fig, [figPath '.png']);
    close(fig);
end

function bar3A_comparison(data, configLabels, designLabels, colors, yLab, titleStr)
    % data: [nConf × nDes]
    nConf = size(data, 1);
    nDes  = size(data, 2);
    hold on;
    bw = 0.8 / nDes;
    offsets = linspace(-0.4 + bw/2, 0.4 - bw/2, nDes);
    for di = 1:nDes
        x = (1:nConf) + offsets(di);
        b = bar(x, data(:, di), bw, 'FaceColor', colors{di}, ...
            'EdgeColor', 'none', 'DisplayName', designLabels{di});
    end
    set(gca, 'XTick', 1:nConf, 'XTickLabel', configLabels, ...
        'XTickLabelRotation', 20);
    legend('Location', 'northwest');
    ylabel(yLab);
    title(titleStr);
    grid on;
end

% stage3A_fullArmVerification.m
% Stage 3A: Full-arm verification of reference-module SIMP designs.
%
% Loads two optimised reference-module topologies (bending-only and
% multi-load shell), expands each to the full Arm-Z arm via segmentToArm,
% then runs frameBasedSolver for five representative configurations.
%
% Designs:
%   A  My_only (vf=0.40, Rfilter ~= 3 FE sizes)
%   B  My+Mz+Ms+Ty+Tz (vf=0.40, Rfilter ~= 3 FE sizes, p=4)
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
arm = armModelDefaults("thin");
E = arm.E;
nu = arm.nu;
R = arm.R;
r = arm.r;
h_seg = arm.h_seg;
alpha = arm.alpha;
res_th = arm.res_th;
Pz = arm.Pz;
ShapeFn = arm.ShapeFn;
useParallel = license('test', 'Distrib_Computing_Toolbox');

%% ---- Load module designs -------------------------------------------------
sweepRoot = fullfile(scriptDir, 'results', 'referenceModuleSIMP_sweep');

fprintf('Loading reference-module designs from Stage 2C sweep...\n');

% Design A: bending-only (My only, vf=0.40, Rfilter ~= 3 FE sizes)
dA = load(fullfile(sweepRoot, 'My_only_vf040_rmin300', 'result.mat'), 'rho_opt');
rhoA = dA.rho_opt;

% Design B: multi-load shell (My+Mz+Ms+Ty+Tz, vf=0.40, Rfilter ~= 3 FE sizes, p=4)
dB = load(fullfile(sweepRoot, ...
    'My_Mz_Ms_Ty_Tz_vf040_rmin300_p04', 'result.mat'), 'rho_opt');
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
fprintf('Parallel verification sweep: %d\n', useParallel);

%% ---- Main sweep ---------------------------------------------------------
% Preallocate metric arrays: [nConf × nDes]
tipDisp_mm  = nan(nConf, nDes);
maxDisp_mm  = nan(nConf, nDes);
maxHM_MPa   = nan(nConf, nDes);
volFracs    = nan(nConf, nDes);
statuses    = strings(nConf, nDes);

jobs = {};
for ci = 1:nConf
    for di = 1:nDes
        jobs{end+1, 1} = struct('ci', ci, 'di', di, 'cfg', configs{ci}, 'des', designs{di}); %#ok<SAGROW>
    end
end
jobResults = cell(numel(jobs), 1);

if useParallel
    parfor ji = 1:numel(jobs)
        job = jobs{ji};
        jobResults{ji} = runStage3AVerificationCase(job.ci, job.di, job.cfg, job.des, ...
            E, nu, h_seg, R, r, res_th, alpha, ShapeFn, Pz, ...
            arm.constEndRing, arm.constMiddleRing, resultRoot);
    end
else
    for ji = 1:numel(jobs)
        job = jobs{ji};
        jobResults{ji} = runStage3AVerificationCase(job.ci, job.di, job.cfg, job.des, ...
            E, nu, h_seg, R, r, res_th, alpha, ShapeFn, Pz, ...
            arm.constEndRing, arm.constMiddleRing, resultRoot);
    end
end

for ji = 1:numel(jobResults)
    jr = jobResults{ji};
    ci = jr.ci;
    di = jr.di;
    tipDisp_mm(ci, di) = jr.tipDisp_mm;
    maxDisp_mm(ci, di) = jr.maxDisp_mm;
    maxHM_MPa(ci, di) = jr.maxHM_MPa;
    volFracs(ci, di) = jr.volFrac;
    statuses(ci, di) = jr.status;
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
savefig(fig1, fullfile(resultRoot, 'comparison_tipDisp.fig'));
close(fig1);

% Figure 2: max HM stress comparison
fig2 = figure('Visible', 'off', 'Name', 'maxHM comparison', ...
    'Position', [100 100 900 420]);
bar3A_comparison(maxHM_MPa, configLabels, designLabels, colors, ...
    'Max HM stress (MPa)', 'Stage 3A: Max Huber-Mises stress by design and configuration');
saveas(fig2, fullfile(resultRoot, 'comparison_maxHM.png'));
savefig(fig2, fullfile(resultRoot, 'comparison_maxHM.fig'));
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
savefig(fig3, fullfile(resultRoot, 'stress_ratio_A_over_B.fig'));
close(fig3);

fprintf('\nComparison figures saved to %s\n', resultRoot);

% =========================================================================
function result = runStage3AVerificationCase(ci, di, cfg, des, ...
        E, nu, h_seg, R, r, res_th, alpha, ShapeFn, Pz, ...
        constEndRing, constMiddleRing, resultRoot)
    fprintf('\n=== Config %d, design %d: %s / %s ===\n', ci, di, cfg.label, des.label);
    result = struct('ci', ci, 'di', di, 'tipDisp_mm', NaN, 'maxDisp_mm', NaN, ...
        'maxHM_MPa', NaN, 'volFrac', NaN, 'status', "failed");
    tic;
    try
        model = ManipulatorModel3D(E, nu, h_seg, R, r, 15, res_th, alpha, ...
            cfg.betas, ShapeFn, true, Pz, constEndRing, constMiddleRing);
        nArmElems = model.analysis.getTotalElemsNumber();

        if isempty(des.rhoSeg)
            x_arm = ones(nArmElems, 1);
        else
            x_arm = model.segmentToArm(des.rhoSeg);
            assert(numel(x_arm) == nArmElems, ...
                'segmentToArm returned %d elements, arm model has %d.', ...
                numel(x_arm), nArmElems);
        end

        [~, ~] = model.frameBasedSolver(x_arm);

        qn = model.analysis.qnodal;
        dispNorm_m = sqrt(sum(qn.^2, 2));
        result.maxDisp_mm = max(dispNorm_m) * 1e3;

        tipIds = find(model.loadSurfaceNodes);
        result.tipDisp_mm = norm(mean(qn(tipIds, :), 1)) * 1e3;

        activeElems = x_arm > 0.5;
        if ~any(activeElems)
            activeElems = true(nArmElems, 1);
        end
        hmGP = model.fe.results.gp.all(13, activeElems, :);
        result.maxHM_MPa = max(hmGP(:)) / 1e6;
        result.volFrac = mean(x_arm);
        result.status = "ok";

        elapsed = toc;
        fprintf('  [%s / %s] tipDisp=%.3f mm  maxHM=%.4f MPa  vf=%.3f  (%.1f s)\n', ...
            cfg.label, des.label, result.tipDisp_mm, result.maxHM_MPa, result.volFrac, elapsed);

        figPath = fullfile(resultRoot, sprintf('%s__%s', cfg.name, des.name));
        saveStressFigure(model, x_arm, activeElems, cfg.label, des.label, ...
            result.maxHM_MPa, result.tipDisp_mm, figPath);
    catch ME
        result.status = "failed";
        fprintf('  [%s / %s] FAILED: %s\n', cfg.label, des.label, ME.message);
        save(fullfile(resultRoot, ...
            sprintf('failed_%s_%s.mat', cfg.name, des.name)), 'ME', 'cfg', 'des');
    end
end

function saveStressFigure(model, x_arm, activeElems, cfgLabel, desLabel, ...
        maxHM_MPa, tipDisp_mm, figPath)
    fig = figure('Visible', 'off', 'Name', [cfgLabel ' - ' desLabel]);
    hold on; axis on; daspect([1 1 1]); view(45, 25);

    % Node-averaged HM stress (MPa) from element results
    hmNodal = model.fe.results.nodal.all(:, 13) / 1e6;

    elems = model.mesh.elems;
    nodes = model.mesh.nodes;
    faceP = model.fe.shapeFn.fcontours';      % [nFacesPerElem × nodesPerFace]
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
    savefig(fig, [figPath '.fig']);
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

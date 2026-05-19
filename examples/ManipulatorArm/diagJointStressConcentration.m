% diagJointStressConcentration
% Diagnostic: are peak HM stresses at segment joints or mid-body?
%
% Four designs are compared, all with the adversarially worst-case beta
% (beta=[0,270,0,90,180,180,180]) and the nominal straight arm (all betas=0):
%
%   1. Full material (x=1), no fillet   - isolates geometric stress concentration
%   2. Helix initial params, no fillet  - approximate current-generation design level
%   3. Nominal beta, x=1, no fillet     - reference loading without large bends
%   4. Full material (x=1), with fillet - quantifies Kt reduction from joint rounding
%
% No optimisation is run; each case is a single FE solve (~10 s).
% Figures are saved to results/diagJointStressConcentration/.

clear; close all; clc;
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

%% ---- Parameters -------------------------------------------------------
arm     = armModelDefaults("thin");
E       = arm.E;
nu      = arm.nu;
R       = arm.R;
r       = arm.r;
h_seg   = arm.h_seg;
alpha   = arm.alpha;
res     = arm.res;
res_th  = arm.res_th;
Pz      = arm.Pz;
ShapeFn = arm.ShapeFn;

% Worst-case beta from CG iteration 1 adversarial sweep (cand 19, s=2.581)
betaWorst   = [0, 270, 0, 90, 180, 180, 180];
% Nominal (straight) arm for reference
betaNominal = [0, 0,   0,  0,   0,   0,   0];

% Fillet radius for Case 4: 30% of wall thickness
wallThickness = R - r;
filletR = 0.3 * wallThickness;
fprintf('  Wall thickness  : %.4f m\n', wallThickness);
fprintf('  Fillet radius   : %.4f m  (%.0f%% of wall)\n\n', filletR, 100*filletR/wallThickness);

resultRoot = fullfile(scriptDir, 'results', 'diagJointStressConcentration');
if ~exist(resultRoot, 'dir'), mkdir(resultRoot); end

fprintf('Joint stress concentration diagnostic\n');
fprintf('  Worst-case beta : %s\n', mat2str(betaWorst));
fprintf('  Nominal beta    : %s\n', mat2str(betaNominal));
fprintf('  Result root     : %s\n\n', resultRoot);

%% ---- Build models -----------------------------------------------------
fprintf('Building models...\n');
modelWorst = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
    betaWorst, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);

modelNominal = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
    betaNominal, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv);

modelFilletWorst = ManipulatorModel3D(E, nu, h_seg, R, r, res, res_th, alpha, ...
    betaWorst, ShapeFn, false, Pz, arm.constEndRing, arm.constMiddleRing, arm.nCircDiv, filletR);

nElems = modelWorst.analysis.getTotalElemsNumber();
fprintf('  Total elements: %d\n\n', nElems);

%% ---- Case 1: full material, worst-case beta ---------------------------
fprintf('Case 1: full material, worst-case beta...\n');
xFull = ones(nElems, 1);
modelWorst.analysis.solveWeighted(xFull);
modelWorst.analysis.computeElementResults(xFull);

hmAll_worst_full = squeeze(modelWorst.analysis.felems{1}.results.gp.all(13, :, :));
hmElem_worst_full = max(hmAll_worst_full, [], 2);   % max over GPs per element
maxHM_worst_full  = max(hmElem_worst_full);
fprintf('  maxHM = %.4e Pa\n', maxHM_worst_full);

fig1 = figure('Name', 'Case 1: full material, worst-case beta', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
light('Position', [1 1 5], 'Style', 'infinite');
modelWorst.analysis.plotMaps("sHM", 0.0);
colorbar; view(3); axis equal off;
title(sprintf('Full material, worst-case \\beta   maxHM=%.3e Pa', maxHM_worst_full));
exportgraphics(fig1, fullfile(resultRoot, 'case1_full_worst_3d.png'), 'Resolution', 200);
savefig(fig1, fullfile(resultRoot, 'case1_full_worst_3d.fig'));

% Front view
fig1b = figure('Name', 'Case 1 front', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
modelWorst.analysis.plotMaps("sHM", 0.0);
colorbar; view(0, 0); axis equal off;
title('Full material, worst-case \beta  (front view)');
exportgraphics(fig1b, fullfile(resultRoot, 'case1_full_worst_front.png'), 'Resolution', 200);
savefig(fig1b, fullfile(resultRoot, 'case1_full_worst_front.fig'));

%% ---- Case 2: helix initial params at VF≈0.51, worst-case beta ---------
fprintf('\nCase 2: helix initial params at VF=0.51, worst-case beta...\n');
curveOpts = struct();
curveOpts.penal          = 3.0;
curveOpts.constRefElems  = armConstRingElementIds(modelWorst, arm, "linked");
curveOpts.constFullElems = armConstRingElementIds(modelWorst, arm, "full");
curveOpts.rhoMin         = arm.mma.xminValue;
curveOpts.elemSize       = arm.nominalElementSize;
curveOpts.VolFrac        = 0.51;

p0 = defaultCurveParamInitial(arm);
[~, xHelix] = buildCurveLinkedDensity(p0, modelWorst, curveOpts);
actualVF = mean(xHelix);
fprintf('  Actual VF = %.4f\n', actualVF);

modelWorst.analysis.solveWeighted(xHelix);
modelWorst.analysis.computeElementResults(xHelix);

hmAll_worst_helix = squeeze(modelWorst.analysis.felems{1}.results.gp.all(13, :, :));
hmElem_worst_helix = max(hmAll_worst_helix, [], 2);
maxHM_worst_helix  = max(hmElem_worst_helix(xHelix > 0.5));
fprintf('  maxHM (x>0.5) = %.4e Pa\n', maxHM_worst_helix);

modelWorst.analysis.felems{1}.selectedElems = xHelix > 0.5;
fig2 = figure('Name', 'Case 2: helix VF=0.51, worst-case beta', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
light('Position', [1 1 5], 'Style', 'infinite');
modelWorst.analysis.plotMaps("sHM", 0.0);
colorbar; view(3); axis equal off;
title(sprintf('Helix VF=%.2f, worst-case \\beta   maxHM=%.3e Pa', actualVF, maxHM_worst_helix));
exportgraphics(fig2, fullfile(resultRoot, 'case2_helix_worst_3d.png'), 'Resolution', 200);
savefig(fig2, fullfile(resultRoot, 'case2_helix_worst_3d.fig'));

fig2b = figure('Name', 'Case 2 front', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
modelWorst.analysis.plotMaps("sHM", 0.0);
colorbar; view(0, 0); axis equal off;
title('Helix VF=0.51, worst-case \beta  (front view)');
exportgraphics(fig2b, fullfile(resultRoot, 'case2_helix_worst_front.png'), 'Resolution', 200);
savefig(fig2b, fullfile(resultRoot, 'case2_helix_worst_front.fig'));

modelWorst.analysis.felems{1}.selectedElems = [];

%% ---- Case 3: full material, nominal beta (reference) ------------------
fprintf('\nCase 3: full material, nominal beta...\n');
modelNominal.analysis.solveWeighted(xFull);
modelNominal.analysis.computeElementResults(xFull);

hmAll_nominal = squeeze(modelNominal.analysis.felems{1}.results.gp.all(13, :, :));
hmElem_nominal = max(hmAll_nominal, [], 2);
maxHM_nominal  = max(hmElem_nominal);
fprintf('  maxHM = %.4e Pa\n', maxHM_nominal);

fig3 = figure('Name', 'Case 3: full material, nominal beta', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
light('Position', [1 1 5], 'Style', 'infinite');
modelNominal.analysis.plotMaps("sHM", 0.0);
colorbar; view(3); axis equal off;
title(sprintf('Full material, nominal \\beta   maxHM=%.3e Pa', maxHM_nominal));
exportgraphics(fig3, fullfile(resultRoot, 'case3_full_nominal_3d.png'), 'Resolution', 200);
savefig(fig3, fullfile(resultRoot, 'case3_full_nominal_3d.fig'));

fig3b = figure('Name', 'Case 3 front', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
modelNominal.analysis.plotMaps("sHM", 0.0);
colorbar; view(0, 0); axis equal off;
title('Full material, nominal \beta  (front view)');
exportgraphics(fig3b, fullfile(resultRoot, 'case3_full_nominal_front.png'), 'Resolution', 200);
savefig(fig3b, fullfile(resultRoot, 'case3_full_nominal_front.fig'));

%% ---- Case 4: full material, worst-case beta, with fillet --------------
fprintf('\nCase 4: full material, worst-case beta, filletR=%.4f m...\n', filletR);
modelFilletWorst.analysis.solveWeighted(xFull);
modelFilletWorst.analysis.computeElementResults(xFull);

hmAll_fillet_full = squeeze(modelFilletWorst.analysis.felems{1}.results.gp.all(13, :, :));
hmElem_fillet_full = max(hmAll_fillet_full, [], 2);
maxHM_fillet_full  = max(hmElem_fillet_full);
fprintf('  maxHM = %.4e Pa  (%.1f%% of no-fillet)\n', maxHM_fillet_full, ...
    100 * maxHM_fillet_full / maxHM_worst_full);

fig4a = figure('Name', 'Case 4: fillet, worst-case beta', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
light('Position', [1 1 5], 'Style', 'infinite');
modelFilletWorst.analysis.plotMaps("sHM", 0.0);
colorbar; view(3); axis equal off;
title(sprintf('Fillet r=%.3fm, worst-case \\beta   maxHM=%.3e Pa', filletR, maxHM_fillet_full));
exportgraphics(fig4a, fullfile(resultRoot, 'case4_fillet_worst_3d.png'), 'Resolution', 200);
savefig(fig4a, fullfile(resultRoot, 'case4_fillet_worst_3d.fig'));

fig4b = figure('Name', 'Case 4 front', 'Position', [100 100 900 600]);
hold on; axis on; daspect([1 1 1]);
light('Position', [-1 -2 5], 'Style', 'local');
modelFilletWorst.analysis.plotMaps("sHM", 0.0);
colorbar; view(0, 0); axis equal off;
title(sprintf('Fillet r=%.3fm, worst-case \\beta  (front view)', filletR));
exportgraphics(fig4b, fullfile(resultRoot, 'case4_fillet_worst_front.png'), 'Resolution', 200);
savefig(fig4b, fullfile(resultRoot, 'case4_fillet_worst_front.fig'));

%% ---- Per-element stress: identify joint vs mid-body peaks -------------
fprintf('\nLocating peak-stress elements...\n');

H          = modelWorst.halfSegmentNelems;
nSegments  = nElems / H;
segmentIdx = ceil((1:nElems)' / H);   % segment index for each element (1-based half-segments)
jointBoundaryFrac = 0.15;             % elements within 15% of segment ends = "joint zone"
localIdx   = mod((1:nElems)' - 1, H) + 1;
inJointZone = localIdx <= round(jointBoundaryFrac * H) | ...
              localIdx >= round((1 - jointBoundaryFrac) * H);

% Case 1 summary
[~, iPeak1] = max(hmElem_worst_full);
fprintf('  Case 1 peak element: %d  segment: %d  in-joint-zone: %d\n', ...
    iPeak1, segmentIdx(iPeak1), inJointZone(iPeak1));

jointMax1   = max(hmElem_worst_full(inJointZone));
midBodyMax1 = max(hmElem_worst_full(~inJointZone));
fprintf('    Joint-zone max  = %.4e Pa  (%.1f%% of global max)\n', ...
    jointMax1,   100 * jointMax1   / maxHM_worst_full);
fprintf('    Mid-body max    = %.4e Pa  (%.1f%% of global max)\n', ...
    midBodyMax1, 100 * midBodyMax1 / maxHM_worst_full);

% Case 3 summary (nominal, for comparison)
[~, iPeak3] = max(hmElem_nominal);
fprintf('  Case 3 peak element: %d  segment: %d  in-joint-zone: %d\n', ...
    iPeak3, segmentIdx(iPeak3), inJointZone(iPeak3));

jointMax3   = max(hmElem_nominal(inJointZone));
midBodyMax3 = max(hmElem_nominal(~inJointZone));
fprintf('    Joint-zone max  = %.4e Pa  (%.1f%% of global max)\n', ...
    jointMax3,   100 * jointMax3   / maxHM_nominal);
fprintf('    Mid-body max    = %.4e Pa  (%.1f%% of global max)\n', ...
    midBodyMax3, 100 * midBodyMax3 / maxHM_nominal);

% Case 4 summary (fillet, worst-case beta)
[~, iPeak4] = max(hmElem_fillet_full);
fprintf('  Case 4 peak element: %d  segment: %d  in-joint-zone: %d\n', ...
    iPeak4, segmentIdx(iPeak4), inJointZone(iPeak4));

jointMax4   = max(hmElem_fillet_full(inJointZone));
midBodyMax4 = max(hmElem_fillet_full(~inJointZone));
fprintf('    Joint-zone max  = %.4e Pa  (%.1f%% of global max)\n', ...
    jointMax4,   100 * jointMax4   / maxHM_fillet_full);
fprintf('    Mid-body max    = %.4e Pa  (%.1f%% of global max)\n', ...
    midBodyMax4, 100 * midBodyMax4 / maxHM_fillet_full);
fprintf('  Fillet Kt reduction: %.3f -> %.3f  (%.1f%% drop in peak)\n', ...
    maxHM_worst_full / midBodyMax1, maxHM_fillet_full / midBodyMax4, ...
    100 * (1 - maxHM_fillet_full / maxHM_worst_full));

%% ---- Summary bar chart ------------------------------------------------
figSum = figure('Name', 'Stress summary', 'Position', [100 100 800 480]);
cases   = {'Full, worst \beta', 'Helix VF=0.51, worst \beta', ...
           'Full, nominal \beta', sprintf('Fillet r=%.3fm, worst \\beta', filletR)};
maxVals = [maxHM_worst_full, maxHM_worst_helix, maxHM_nominal, maxHM_fillet_full] / 1e6;
colors  = [0.28 0.45 0.72; 0.28 0.45 0.72; 0.55 0.72 0.42; 0.85 0.55 0.25];
b = bar(maxVals, 'FaceColor', 'flat');
b.CData = colors;
set(gca, 'XTickLabel', cases, 'XTick', 1:4);
xtickangle(15);
ylabel('Peak HM stress [MPa]');
title('Peak HM stress comparison');
yline(maxHM_worst_full / 1e6, '--k', 'No-fillet ref', 'LabelHorizontalAlignment', 'left');
grid on;
exportgraphics(figSum, fullfile(resultRoot, 'summary_peak_stress.png'), 'Resolution', 200);
savefig(figSum, fullfile(resultRoot, 'summary_peak_stress.fig'));

%% ---- Save CSV summary -------------------------------------------------
T = table( ...
    {'full_worst'; 'helix_worst'; 'full_nominal'; 'fillet_worst'}, ...
    [maxHM_worst_full; maxHM_worst_helix; maxHM_nominal; maxHM_fillet_full], ...
    [jointMax1; NaN; jointMax3; jointMax4], ...
    [midBodyMax1; NaN; midBodyMax3; midBodyMax4], ...
    'VariableNames', {'case', 'maxHM_Pa', 'jointZoneMax_Pa', 'midBodyMax_Pa'});
writetable(T, fullfile(resultRoot, 'stress_summary.csv'));

fprintf('\nDiagnostic complete. Figures saved to:\n  %s\n', resultRoot);
fprintf('\nInterpretation guide:\n');
fprintf('  jointZoneMax / maxHM > 0.90  ->  joint geometry drives the peak\n');
fprintf('  jointZoneMax / maxHM < 0.70  ->  mid-body bending dominates\n');
fprintf('  Case4/Case1 ratio < 0.78     ->  fillet breaks stress floor (Kt<1.1)\n');

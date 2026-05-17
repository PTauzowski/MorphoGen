% diagJointStressConcentration
% Diagnostic: are peak HM stresses at segment joints or mid-body?
%
% Three designs are compared, all with the adversarially worst-case beta
% (cand 19 from the CG iter-1 sweep: beta=[0,270,0,90,180,180,180],
% s_ratio=2.581) and the nominal straight arm (all betas=0):
%
%   1. Full material (x=1)    - isolates geometric stress concentration
%   2. Helix initial params   - approximate current-generation design level
%   3. Nominal beta, x=1      - reference loading without large bends
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

%% ---- Summary bar chart ------------------------------------------------
fig4 = figure('Name', 'Stress summary', 'Position', [100 100 700 450]);
cases   = {'Full, worst \beta', 'Helix VF=0.51, worst \beta', 'Full, nominal \beta'};
maxVals = [maxHM_worst_full, maxHM_worst_helix, maxHM_nominal] / 1e6;
bar(maxVals, 'FaceColor', [0.28 0.45 0.72]);
set(gca, 'XTickLabel', cases, 'XTick', 1:3);
ylabel('Peak HM stress [MPa]');
title('Peak HM stress comparison');
grid on;
exportgraphics(fig4, fullfile(resultRoot, 'summary_peak_stress.png'), 'Resolution', 200);
savefig(fig4, fullfile(resultRoot, 'summary_peak_stress.fig'));

%% ---- Save CSV summary -------------------------------------------------
T = table( ...
    {'full_worst'; 'helix_worst'; 'full_nominal'}, ...
    [maxHM_worst_full; maxHM_worst_helix; maxHM_nominal], ...
    [jointMax1; NaN; jointMax3], ...
    [midBodyMax1; NaN; midBodyMax3], ...
    'VariableNames', {'case', 'maxHM_Pa', 'jointZoneMax_Pa', 'midBodyMax_Pa'});
writetable(T, fullfile(resultRoot, 'stress_summary.csv'));

fprintf('\nDiagnostic complete. Figures saved to:\n  %s\n', resultRoot);
fprintf('\nInterpretation guide:\n');
fprintf('  jointZoneMax / maxHM > 0.90  ->  joint geometry drives the peak\n');
fprintf('  jointZoneMax / maxHM < 0.70  ->  mid-body bending dominates\n');

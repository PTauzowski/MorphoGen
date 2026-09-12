%% coupledExamples.m
% Generate illustrative comparison figures for the coupled frame-solid
% solver and the corresponding full 3D solid analysis. The script exports
% paired deformed-mesh and Huber-Mises stress figures for representative
% joint-angle configurations and reports the frame cross-section constants
% used by the beam surrogate.

baseDir = fileparts(mfilename('fullpath'));
repoDir = fullfile(baseDir, '..', '..');

addpath(genpath(fullfile(repoDir, 'analysis')));
addpath(genpath(fullfile(repoDir, 'design')));
addpath(genpath(fullfile(repoDir, 'elements')));
addpath(genpath(fullfile(repoDir, 'examples', 'models')));
addpath(genpath(fullfile(repoDir, 'materials')));
addpath(genpath(fullfile(repoDir, 'math')));
addpath(genpath(fullfile(repoDir, 'mesh')));
addpath(genpath(fullfile(repoDir, 'examples', 'ManipulatorArm', 'HelperFunctions')));
addpath(genpath(fullfile(repoDir, 'postprocessing')));

docsDir = fullfile(repoDir, 'docs');
if ~exist(docsDir, 'dir')
    mkdir(docsDir);
end

close all;

arm = armModelDefaults("thin");
E = arm.E;
nu = arm.nu;
ShapeFn = arm.ShapeFn;
R = arm.R;
r = arm.r;
alpha = arm.alpha;
segmentLength = arm.segmentLength;
res = arm.res;
res_thickness = arm.res_th;
Pz = arm.Pz;

frameElems=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
nArms = size(frameElems,1);
nDiv  = 8;                           % 8 bins per angle
%vals  = linspace(0, 360 - 360/nDiv, nDiv);  % [0,45,90,...,315]
vals  = linspace(-180, 180 - 360/nDiv, nDiv);  % [0,45,90,...,315]

% We keep the first column fixed at 0, vary the remaining 6
nVar = nArms - 1;                    % 6
G = cell(1, nVar);
[G{:}] = ndgrid(vals);

% Assemble samples: rows = 8^6, cols = 7 (first col fixed to 0)
samples = zeros(nDiv^nVar, nArms);
%samples(:,1) = 0;
for k = 1:nVar
    samples(:, k+1) = G{k}(:);
end

nSamples = size(samples,1);
[~, ~, frameNodes] = computeArmSamples(E, nu, segmentLength, alpha, samples );
[vN, vTy, vTz, vMs, vMy, vMz, all_forces] = computeAllInternalForces( frameNodes, frameElems, E, nu, R, r, samples, Pz);

save("ManipulatorFrame3Ddata200K.mat", "-v7.3");
%load("ManipulatorFrame3Ddata200K.mat");

% --- Statistically-derived extremal configurations -----------------------
% Extract the grid sample that maximises each force quantity.  The mirror
% symmetry  forces(-beta) = -forces(beta)  is preserved: min_* = -max_*.
extremalCfgs = extractExtremalConfigs(vN, vTy, vTz, vMs, vMy, vMz, samples);

% Save for use by armLoadConfigs("statistical") in optimisation scripts.
dataDir = fullfile(baseDir, 'data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end
configs = extremalCfgs;  %#ok<NASGU>
save(fullfile(dataDir, 'extremalLoadConfigs.mat'), 'configs');
fprintf('Saved extremal load configs to %s\n', fullfile(dataDir, 'extremalLoadConfigs.mat'));

% Build visualization examples from the four positive extremals.
% "min_*" configs are just the mirrored arm; showing max_* covers
% all distinct structural situations.
positiveMask = cellfun(@(c) ~startsWith(c.name, 'min'), extremalCfgs);
visCfgs = extremalCfgs(positiveMask);
captionMap = containers.Map( ...
    {'max_bending', 'max_torsion', 'max_shear', 'max_tension'}, ...
    {'Max bending-moment config (extremal M_b)', ...
     'Max torsion config (extremal M_s)', ...
     'Max shear config (extremal T)', ...
     'Max axial-tension config (extremal N)'});
examples = struct( ...
    'name',    cellfun(@(c) c.name,             visCfgs, 'UniformOutput', false), ...
    'sample',  cellfun(@(c) c.betas,            visCfgs, 'UniformOutput', false), ...
    'caption', cellfun(@(c) captionMap(c.name), visCfgs, 'UniformOutput', false) ...
);

nbin=500;

figure;
histogram(vN(:,1,2),nbin);
title('N 1');

figure;
histogram(vTy(:,1,2),nbin);
title('Ty 1');

figure;
histogram(vTz(:,1,2),nbin);
title('Tz 1');

figure;
histogram(vMs(:,1,2),nbin);
title('Ms 1');

figure;
histogram(vMy(:,1,2),nbin);
title('My 1');

figure;
histogram(vMz(:,1,2),nbin);
title('Mz 1');

figure;
histogram(vN(:,2,2),nbin);
title('N 2');

figure;
histogram(vTy(:,2,2),nbin);
title('Ty 2');

figure;
histogram(vTz(:,2,2),nbin);
title('Tz 2');

figure;
histogram(vMs(:,2,2),nbin);
title('Ms 2');

figure;
histogram(vMy(:,2,2),nbin);
title('My 2');

figure;
histogram(vMz(:,2,2),nbin);
title('Mz 2');

Ms_env = max(max(abs(vMs), [], 2), [], 3);
T_env  = max(max(hypot(vTy, vTz), [], 2), [], 3);
Mb_env = max(max(hypot(vMy, vMz), [], 2), [], 3);

figure; histogram(Ms_env, nbin, 'Normalization', 'probability'); title('|M_s| envelope');
figure; histogram(T_env, nbin, 'Normalization', 'probability'); title('shear resultant envelope');
figure; histogram(Mb_env, nbin, 'Normalization', 'probability'); title('bending resultant envelope');

sectionReportPrinted = false;

for k = 1:numel(examples)
    example = examples(k);
    fprintf('Generating comparison example: %s\n', example.name);

    coupledModel = ManipulatorModel3D( ...
        E, nu, segmentLength, R, r, res, res_thickness, alpha, ...
        example.sample, ShapeFn, true, Pz, arm.constEndRing, arm.constMiddleRing);
    [qCoupled, ~] = coupledModel.frameBasedSolver(1);

    directModel = ManipulatorModel3D( ...
        E, nu, segmentLength, R, r, res, res_thickness, alpha, ...
        example.sample, ShapeFn, true, Pz, arm.constEndRing, arm.constMiddleRing);
    directModel.compute(1);
    qDirect = directModel.analysis.qnodal;

    if ~sectionReportPrinted
        reportFrameSectionConsistency(coupledModel, E, nu, R, r);
        sectionReportPrinted = true;
    end

    exportDeformedComparisonPlot( ...
        coupledModel, qCoupled, directModel, qDirect, example, docsDir);
    exportStressComparisonPlot(coupledModel, directModel, example, docsDir);
    exportStressComparisonNormalized(coupledModel, directModel, example, docsDir);
    close all;
end

function exportDeformedComparisonPlot(coupledModel, qCoupled, directModel, qDirect, example, docsDir)
    fig = figure('Color', 'w', 'Position', [100 100 1600 720]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    maxDispCoupled = max(vecnorm(qCoupled, 2, 2));
    maxDispDirect = max(vecnorm(qDirect, 2, 2));
    globalMaxDisp = max([maxDispCoupled; maxDispDirect; eps]);

    ax1 = nexttile(tl, 1);
    defNodesCoupled = plotDeformedPanel(ax1, coupledModel, qCoupled, globalMaxDisp, ...
        'Coupled frame-solid', maxDispCoupled);

    ax2 = nexttile(tl, 2);
    defNodesDirect = plotDeformedPanel(ax2, directModel, qDirect, globalMaxDisp, ...
        'Full 3D solid', maxDispDirect);

    syncAxes([ax1 ax2], [coupledModel.mesh.nodes; directModel.mesh.nodes; defNodesCoupled; defNodesDirect]);
    sgtitle(strrep(example.name, '_', '\_'), 'Interpreter', 'tex');

    exportgraphics( ...
        fig, fullfile(docsDir, ['coupledExample_' example.name '_deformed.png']), ...
        'Resolution', 600);
    close(fig);
end

function exportStressComparisonPlot(coupledModel, directModel, example, docsDir)
    % Shared color scale — both panels use [globalMin, globalMax].
    % This reveals absolute stress magnitude differences but may suppress
    % distribution patterns in the lower-stress panel.
    fig = figure('Color', 'w', 'Position', [100 100 1600 720]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    hmCoupled = coupledModel.fe.results.nodal.all(:, 13);
    hmDirect  = directModel.fe.results.nodal.all(:, 13);
    hmLimits  = [min([hmCoupled; hmDirect]), max([hmCoupled; hmDirect])];

    ax1 = nexttile(tl, 1);
    defNodesCoupled = plotStressPanel(ax1, coupledModel, 'Coupled frame-solid', hmLimits);

    ax2 = nexttile(tl, 2);
    defNodesDirect = plotStressPanel(ax2, directModel, 'Full 3D solid', hmLimits);

    syncAxes([ax1 ax2], [defNodesCoupled; defNodesDirect]);
    sgtitle([strrep(example.name, '_', '\_') ' — shared color scale'], 'Interpreter', 'tex');

    exportgraphics( ...
        fig, fullfile(docsDir, ['coupledExample_' example.name '_stress.png']), ...
        'Resolution', 600);
    close(fig);
end

function exportStressComparisonNormalized(coupledModel, directModel, example, docsDir)
    % Per-panel (independent) color scale — each panel is normalized to its
    % own [0, HM_max].  This reveals whether the spatial stress distribution
    % patterns are qualitatively similar even when the magnitudes differ.
    fig = figure('Color', 'w', 'Position', [100 100 1600 720]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    maxHM_coupled = max(reshape(coupledModel.fe.results.gp.all(13, :, :), [], 1));
    maxHM_direct  = max(reshape(directModel.fe.results.gp.all(13, :, :),  [], 1));
    limCoupled = [0, maxHM_coupled];
    limDirect  = [0, maxHM_direct];

    ax1 = nexttile(tl, 1);
    defNodesCoupled = plotStressPanel(ax1, coupledModel, 'Coupled frame-solid', limCoupled);

    ax2 = nexttile(tl, 2);
    defNodesDirect = plotStressPanel(ax2, directModel, 'Full 3D solid', limDirect);

    syncAxes([ax1 ax2], [defNodesCoupled; defNodesDirect]);
    sgtitle([strrep(example.name, '_', '\_') ' — independent color scales'], 'Interpreter', 'tex');

    exportgraphics( ...
        fig, fullfile(docsDir, ['coupledExample_' example.name '_stress_normalized.png']), ...
        'Resolution', 600);
    close(fig);
end

function defNodes = plotDeformedPanel(ax, model, qnodal, globalMaxDisp, panelLabel, maxDisp)
    axes(ax);
    hold(ax, 'on');
    axis(ax, 'on');
    daspect(ax, [1 1 1]);
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    zlabel(ax, 'z');

    plotFaces = getExteriorPlotFaces(model.fe);
    baseNodes = model.mesh.nodes;
    dg = norm(max(baseNodes) - min(baseNodes));
    deformScale = 0.35;

    if globalMaxDisp <= eps
        defNodes = baseNodes;
    else
        defNodes = baseNodes + qnodal / globalMaxDisp * dg * deformScale;
    end

    patch(ax, 'Vertices', baseNodes, 'Faces', plotFaces, ...
        'FaceColor', [0.82 0.82 0.82], 'EdgeColor', 'none', 'FaceAlpha', 0.12);
    patch(ax, 'Vertices', defNodes, 'Faces', plotFaces, ...
        'FaceColor', [0.80 0.16 0.16], 'EdgeColor', 'none', 'FaceAlpha', 0.92);

    light(ax, 'Position', [-1 -2 5], 'Style', 'local');
    light(ax, 'Position', [1 1 5], 'Style', 'infinite');
    lighting(ax, 'gouraud');
    view(ax, 45, 25);
    ax.FontSize = 16;

    title(ax, { ...
        panelLabel, ...
        sprintf('Deformed mesh, u_{max} = %.2f mm', maxDisp * 1000) ...
    }, 'Interpreter', 'tex');
end

function defNodes = plotStressPanel(ax, model, panelLabel, hmLimits)
    axes(ax);
    hold(ax, 'on');

    model.fe.plotMap(model.mesh.nodes, model.analysis.qnodal, 13, 0.2);
    axis(ax, 'on');
    daspect(ax, [1 1 1]);
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    zlabel(ax, 'z');
    caxis(ax, hmLimits);
    view(ax, 45, 25);
    ax.FontSize = 16;

    maxHM = max(reshape(model.fe.results.gp.all(13, :, :), [], 1));
    title(ax, { ...
        panelLabel, ...
        sprintf('Huber-Mises stress, HM_{max} = %.2f kPa', maxHM / 1000) ...
    }, 'Interpreter', 'tex');

    defNodes = model.mesh.nodes + 0.2 * model.analysis.qnodal;
end

function plotFaces = getExteriorPlotFaces(fe)
    allFaces = reshape(fe.elems(:, fe.sf.fcontours)', ...
        size(fe.sf.fcontours, 1), size(fe.sf.fcontours, 2) * size(fe.elems, 1))';
    [~, ifaces] = unique(sort(allFaces, 2), 'rows');
    outerFaces = allFaces(ifaces, :);
    duplicateFaces = allFaces;
    duplicateFaces(ifaces, :) = [];
    [~, exteriorIdx, ~] = setxor(sort(outerFaces, 2), sort(duplicateFaces, 2), 'rows');
    plotFaces = outerFaces(exteriorIdx, :);
end

function syncAxes(axs, nodes)
    mins = min(nodes, [], 1);
    maxs = max(nodes, [], 1);
    span = max(maxs - mins);
    if span <= eps
        span = 1;
    end
    margin = 0.05 * span;

    for k = 1:numel(axs)
        xlim(axs(k), [mins(1) - margin, maxs(1) + margin]);
        ylim(axs(k), [mins(2) - margin, maxs(2) + margin]);
        zlim(axs(k), [mins(3) - margin, maxs(3) + margin]);
    end
end

function reportFrameSectionConsistency(modelRef, E, nu, R, r)
    frameElem = modelRef.frameElem;

    expected.A = pi * (R^2 - r^2);
    expected.G = E / (2 * (1 + nu));
    expected.Iy = (pi / 4) * (R^4 - r^4);
    expected.Iz = expected.Iy;
    expected.Jt = (pi / 2) * (R^4 - r^4);

    relErr = @(actual, ref) abs(actual - ref) / max(abs(ref), eps);

    fprintf('\nFrame section consistency check\n');
    fprintf('  E   = %.6e Pa\n', frameElem.E);
    fprintf('  G   = %.6e Pa (expected %.6e, rel.err %.3e)\n', ...
        frameElem.G, expected.G, relErr(frameElem.G, expected.G));
    fprintf('  A   = %.6e m^2 (expected %.6e, rel.err %.3e)\n', ...
        frameElem.A, expected.A, relErr(frameElem.A, expected.A));
    fprintf('  Iy  = %.6e m^4 (expected %.6e, rel.err %.3e)\n', ...
        frameElem.Jy, expected.Iy, relErr(frameElem.Jy, expected.Iy));
    fprintf('  Iz  = %.6e m^4 (expected %.6e, rel.err %.3e)\n', ...
        frameElem.Jz, expected.Iz, relErr(frameElem.Jz, expected.Iz));
    fprintf('  Jt  = %.6e m^4 (expected %.6e, rel.err %.3e)\n', ...
        frameElem.Ks, expected.Jt, relErr(frameElem.Ks, expected.Jt));
    fprintf('  nu  = %.4f\n', nu);
    fprintf('  note: static frame solve uses E, G, A, Iy, Iz, Jt; density is not used here.\n\n');

    assert(relErr(frameElem.G, expected.G) < 1.0e-12, ...
        'Frame shear modulus G is inconsistent with the solid material.');
    assert(relErr(frameElem.A, expected.A) < 1.0e-12, ...
        'Frame area A is inconsistent with the hollow solid section.');
    assert(relErr(frameElem.Jy, expected.Iy) < 1.0e-12, ...
        'Frame Iy is inconsistent with the hollow solid section.');
    assert(relErr(frameElem.Jz, expected.Iz) < 1.0e-12, ...
        'Frame Iz is inconsistent with the hollow solid section.');
    assert(relErr(frameElem.Ks, expected.Jt) < 1.0e-12, ...
        'Frame torsion constant is inconsistent with the hollow solid section.');
end

% stage3A_extendedMetrics.m
% Extended metrics for Stage 3A full-arm verification.
%
% For each (config × design) pair this script:
%   1. Loads an existing per-case .mat file if it contains the required
%      solution fields (x_arm, elemHM_Pa, elemSED, qnodal_m).
%   2. Otherwise rebuilds the ManipulatorModel3D, calls frameBasedSolver,
%      and saves a compact .mat file for future reuse.
%
% Metrics computed at two active-element thresholds (rho > 0.5 and rho > 0.3):
%   - max HM stress (MPa)
%   - 99th / 95th percentile HM stress (MPa)
%   - mean HM stress on active elements (MPa)
%   - total strain energy (J) and SED fraction on active elements
%   - active element count and active volume fraction
%   - max HM normalised by mean density
%
% Outputs:
%   results/stage3A_fullArmVerification/extended_summary.csv
%   results/stage3A_fullArmVerification/ext_maxHM_*.png  (bar charts)
%   results/stage3A_fullArmVerification/ext_pct99HM_*.png
%   per-case  <cfg.name>__<des.name>.mat  (saved on first solve)

clear; close all; clc;
clear classes;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..', '..');
addpath(genpath(projectRoot));

rng(0, 'twister');

%% ---- Geometry (must match Stage 2 / stage3A) ----------------------------
E      = 2.0e9;
nu     = 0.35;
R      = 0.14;
r      = 0.08;
h_seg  = 0.25;
alpha  = 22.5;
res_th = 4;
Pz     = 100;
ShapeFn = ShapeFunctionL8();

%% ---- Load module designs ------------------------------------------------
sweepRoot = fullfile(scriptDir, 'results', 'referenceModuleSIMP_sweep');

dA   = load(fullfile(sweepRoot, 'My_only_vf040_rmin050', 'result.mat'), 'rho_opt');
rhoA = dA.rho_opt;
dB   = load(fullfile(sweepRoot, ...
    'My_Mz_Ms_Ty_Tz_vf040_rmin050_p04_vf040_rmin050_p04', 'result.mat'), 'rho_opt');
rhoB = dB.rho_opt;

designs = {
    struct('name', 'A_My_only',        'label', 'My only',         'rhoSeg', rhoA, 'color', [0.20 0.55 1.00]);
    struct('name', 'B_multiload_shell', 'label', 'Multi-load shell','rhoSeg', rhoB, 'color', [1.00 0.40 0.15]);
    struct('name', 'C_full_solid',      'label', 'Full solid',      'rhoSeg', [],   'color', [0.35 0.75 0.35]);
};
nDes = numel(designs);

%% ---- Arm configurations -------------------------------------------------
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

%% ---- Required fields to consider a .mat cache valid --------------------
REQUIRED_FIELDS = {'x_arm', 'elemHM_Pa', 'elemSED', 'qnodal_m'};

%% ---- Active-element thresholds ------------------------------------------
thresholds = [0.5, 0.3];
threshNames = {'gt05', 'gt03'};

%% ---- Preallocate result tables ------------------------------------------
% One struct per (config × design × threshold) row
rows = {};

%% ---- Main sweep ---------------------------------------------------------
for ci = 1:nConf
    cfg = configs{ci};
    fprintf('\n=== Config %d/%d: %s ===\n', ci, nConf, cfg.label);

    for di = 1:nDes
        des = designs{di};
        caseTag = sprintf('%s__%s', cfg.name, des.name);
        matFile = fullfile(resultRoot, [caseTag '.mat']);

        %% -- Load or solve ------------------------------------------------
        needSolve = true;
        if exist(matFile, 'file')
            tmp = load(matFile);
            if all(isfield(tmp, REQUIRED_FIELDS))
                x_arm     = tmp.x_arm;
                elemHM_Pa = tmp.elemHM_Pa;
                elemSED   = tmp.elemSED;
                qnodal_m  = tmp.qnodal_m;
                needSolve = false;
                fprintf('  [%s] loaded from cache.\n', des.label);
            end
        end

        if needSolve
            fprintf('  [%s] solving ... ', des.label);
            tic;
            try
                model = ManipulatorModel3D(E, nu, h_seg, R, r, 15, res_th, alpha, ...
                    cfg.betas, ShapeFn, true, Pz);
                nArmElems = model.analysis.getTotalElemsNumber();

                if isempty(des.rhoSeg)
                    x_arm = ones(nArmElems, 1);
                else
                    x_arm = model.segmentToArm(des.rhoSeg);
                end

                [~, ~] = model.frameBasedSolver(x_arm);

                %% Extract element HM (mean over GPs): [nElems × 1]
                gpHM = squeeze(model.fe.results.gp.all(13, :, :));  % [nElems × nGP]
                elemHM_Pa = mean(gpHM, 2);                          % [nElems × 1]

                %% Strain energy density per element (mean over GPs): [nElems × 1]
                gpStress = model.fe.results.gp.stress;  % [nElems × nGP × 6]
                gpStrain = model.fe.results.gp.strain;  % [nElems × nGP × 6]
                gpSED    = 0.5 * sum(gpStress .* gpStrain, 3);  % [nElems × nGP]
                elemSED  = mean(gpSED, 2);                       % [nElems × 1]

                %% Nodal displacements
                qnodal_m = model.analysis.qnodal;  % [nNodes × 3]

                elapsed = toc;
                fprintf('done (%.1f s)\n', elapsed);

                %% Save cache
                save(matFile, 'x_arm', 'elemHM_Pa', 'elemSED', 'qnodal_m', ...
                    'cfg', 'des', '-v7.3');
                fprintf('  [%s] cache saved to %s\n', des.label, matFile);

            catch ME
                fprintf('FAILED: %s\n', ME.message);
                save(fullfile(resultRoot, sprintf('failed_%s_%s.mat', cfg.name, des.name)), ...
                    'ME', 'cfg', 'des');
                continue;
            end
        end

        %% -- Compute metrics for each threshold ---------------------------
        nElems   = numel(x_arm);
        totalSED = sum(elemSED);

        for ti = 1:numel(thresholds)
            thr  = thresholds(ti);
            mask = x_arm > thr;
            if ~any(mask), mask = true(nElems, 1); end

            hmActive  = elemHM_Pa(mask) / 1e6;    % MPa
            sedActive = elemSED(mask);

            row.configName  = cfg.name;
            row.configLabel = cfg.label;
            row.designName  = des.name;
            row.designLabel = des.label;
            row.threshold   = thr;
            row.threshName  = threshNames{ti};

            row.volFrac        = mean(x_arm);
            row.activeCount    = sum(mask);
            row.activeVolFrac  = sum(x_arm(mask)) / nElems;

            row.maxHM_MPa     = max(hmActive);
            row.pct99HM_MPa   = prctile(hmActive, 99);
            row.pct95HM_MPa   = prctile(hmActive, 95);
            row.meanHM_MPa    = mean(hmActive);

            row.totalSED_J    = totalSED;
            row.activeSED_J   = sum(sedActive);
            row.activeSEDfrac = sum(sedActive) / (totalSED + eps);

            % Normalise peak HM by volume fraction (same material budget)
            row.maxHM_normVF  = row.maxHM_MPa / (row.volFrac + eps);
            row.pct99HM_normVF = row.pct99HM_MPa / (row.volFrac + eps);

            % Tip displacement
            dispNorm = sqrt(sum(qnodal_m.^2, 2));
            row.maxDisp_mm = max(dispNorm) * 1e3;

            rows{end+1} = row; %#ok<SAGROW>
        end
    end
end

%% ---- Write extended summary CSV -----------------------------------------
T = struct2table(vertcat(rows{:}));
outCSV = fullfile(resultRoot, 'extended_summary.csv');
writetable(T, outCSV);
fprintf('\nExtended summary saved to %s\n', outCSV);

%% ---- Console table (gt05 threshold only) --------------------------------
mask05 = strcmp(T.threshName, 'gt05');
T05 = T(mask05, :);
fprintf('\n--- Extended metrics (active rho > 0.5) ---\n');
disp(T05(:, {'configLabel','designLabel','volFrac','activeCount','maxHM_MPa', ...
    'pct99HM_MPa','meanHM_MPa','activeSEDfrac','maxHM_normVF'}));

%% ---- Bar-chart comparison figures ---------------------------------------
configLabels = cellfun(@(c) c.label, configs, 'UniformOutput', false);
designLabels = cellfun(@(d) d.label, designs, 'UniformOutput', false);
colors       = cellfun(@(d) d.color,  designs, 'UniformOutput', false);

metricSpecs = {
    'maxHM_MPa',      'Max HM stress (MPa)',        'max';
    'pct99HM_MPa',    '99th pct HM stress (MPa)',   'p99';
    'pct95HM_MPa',    '95th pct HM stress (MPa)',   'p95';
    'meanHM_MPa',     'Mean HM stress (MPa)',        'mean';
    'maxHM_normVF',   'Max HM / vol-frac (MPa)',     'maxNormVF';
    'activeSEDfrac',  'Active SED fraction',         'SEDfrac';
};

for thi = 1:numel(thresholds)
    thrName = threshNames{thi};
    maskThr = strcmp(T.threshName, thrName);
    Tthr    = T(maskThr, :);

    for mi = 1:size(metricSpecs, 1)
        col      = metricSpecs{mi, 1};
        ylab     = metricSpecs{mi, 2};
        shortTag = metricSpecs{mi, 3};

        % Reshape to [nConf × nDes]
        data = nan(nConf, nDes);
        for ci = 1:nConf
            for di = 1:nDes
                rowMask = strcmp(Tthr.configName, configs{ci}.name) & ...
                          strcmp(Tthr.designName, designs{di}.name);
                if any(rowMask)
                    data(ci, di) = Tthr.(col)(rowMask);
                end
            end
        end

        fig = figure('Visible', 'off', 'Position', [100 100 920 430]);
        extBarComparison(data, configLabels, designLabels, colors, ylab, ...
            sprintf('Stage 3A [rho>%g]: %s', thresholds(thi), ylab));
        fname = fullfile(resultRoot, ...
            sprintf('ext_%s_%s.png', shortTag, thrName));
        saveas(fig, fname);
        close(fig);
    end
end

fprintf('\nComparison figures saved to %s\n', resultRoot);

% =========================================================================
function extBarComparison(data, configLabels, designLabels, colors, yLab, titleStr)
    nConf = size(data, 1);
    nDes  = size(data, 2);
    hold on;
    bw      = 0.8 / nDes;
    offsets = linspace(-0.4 + bw/2, 0.4 - bw/2, nDes);
    for di = 1:nDes
        x = (1:nConf) + offsets(di);
        bar(x, data(:, di), bw, 'FaceColor', colors{di}, ...
            'EdgeColor', 'none', 'DisplayName', designLabels{di});
    end
    set(gca, 'XTick', 1:nConf, 'XTickLabel', configLabels, ...
        'XTickLabelRotation', 20);
    legend('Location', 'northwest');
    ylabel(yLab);
    title(titleStr);
    grid on;
end

function postResult = postprocessSIMPResult(optResult, model, analyses, Wfilter, opts)
% POSTPROCESSSIMPRESULT  Post-process a SIMP compliance optimisation result.
%
%   postResult = postprocessSIMPResult(optResult, model, analyses, Wfilter, opts)
%
%   Applies density-threshold extraction strategies to the optimised density
%   field and compares them by continuous diagnostic objective or binary FE
%   reanalysis (expensive, opt-in).  Saves plots and a summary CSV if
%   opts.resultRoot is provided.
%
%   Inputs:
%     optResult  struct from solveSIMPComplianceVolumeMMA or
%                solveSIMPVolumeComplianceMMA.  Required fields:
%                  .zFinal   [H x 1]       reference half-segment design variable,
%                                         or full-arm density for unlinked runs
%                  .xFinal   [nElems x 1]  full-arm density
%                  .C0       [nConfigs x 1]  compliance normalisation
%     model      ManipulatorModel3D  (segmentToArm, mesh)
%     analyses   {nConfigs x 1}  FEAnalysis objects
%     Wfilter    [H x H]  filter matrix used during optimisation
%     opts       struct  (all fields optional)
%
%   opts fields:
%     penal        SIMP exponent                      (required — no default)
%     pAgg         p-norm aggregation exponent         (required — no default)
%     weights      [nConfigs x 1] config weights       (default: uniform)
%     VolFrac      target volume fraction for
%                  volume-preserving threshold          (default: mean(xFinal))
%     skipHeaviside  logical — skip Heaviside sharpening   (default: false)
%                  Set true for sensitivity-filtered SIMP (solveSIMPComplianceVolumeMMA)
%                  where z IS the physical density; Wfilter*z has no meaning
%                  for density sharpening.  Use false only for projection-based
%                  SIMP (solveSIMPVolumeStressMMA with useProjection=true).
%     betaSeq      Heaviside beta values to evaluate   (default: [8 16 32 64])
%     eta          Heaviside projection threshold       (default: 0.5)
%     fixedVars    indices forced to 1 in rho-space    (default: [])
%                  (const ring elements — prevents filter blurring at boundaries)
%     runReanalysis  logical — run binary FE reanalysis (default: false)
%                    ~5 min per method for the full arm; use after optimisation
%                    is finalised, not during development
%     runSweep     logical — run threshold sweep with
%                  binary FE reanalysis                (default: false)
%     thresholds   [nT x 1] for sweep                  (default: linspace(0.1,0.9,9))
%     useParallel  logical                              (default: false)
%     resultRoot   path for saving outputs              (default: '' = no saving)
%     configNames  {nConfigs x 1} string labels for CSV
%
%   Output: postResult struct with fields
%     volumeThreshold   struct  volume-preserving extraction
%     heaviside(i)      struct  Heaviside at betaSeq(i)
%     best              struct  selected candidate. Uses J_binary when
%                       runReanalysis or runSweep is true, otherwise uses
%                       J_continuous as a diagnostic fallback.
%     sweep             struct  from sweepThresholdsReanalysis (if runSweep)
%
%   Each candidate struct has:
%     .xContinuous  [nElems x 1] double   continuous density used diagnostically
%     .solidRaw     [nElems x 1] logical  threshold mask before cleanup
%     .solid        [nElems x 1] logical  binary mask after cleanup
%     .xBinary      [nElems x 1] double   solid + xVoid*(~solid)
%     .volFracRaw   scalar
%     .volFrac      scalar
%     .J_continuous scalar  SIMP objective at xContinuous
%     .J_binary     scalar  binary FE objective (only if reanalysis is enabled)
%     .selectedScore scalar
%     .selectedBy   string
%     .label         string

    %% -- Defaults --------------------------------------------------------
    assert(isfield(opts, 'penal'), 'postprocessSIMPResult: opts.penal is required.');
    assert(isfield(opts, 'pAgg'),  'postprocessSIMPResult: opts.pAgg is required.');

    nConfigs    = numel(analyses);
    penal       = opts.penal;
    pAgg        = opts.pAgg;
    weights     = optField(opts, 'weights',       ones(nConfigs,1)/nConfigs);
    VolFrac     = optField(opts, 'VolFrac',        mean(optResult.xFinal));
    betaSeq     = optField(opts, 'betaSeq',        [8 16 32 64]);
    eta         = optField(opts, 'eta',            0.5);
    fixedVars   = optField(opts, 'fixedVars',      []);
    skipHeaviside = optField(opts, 'skipHeaviside', false);
    runReanalysis = optField(opts, 'runReanalysis', false);
    runSweep    = optField(opts, 'runSweep',        false);
    thresholds  = optField(opts, 'thresholds',      linspace(0.1, 0.9, 9)');
    useParallel = optField(opts, 'useParallel',     false);
    resultRoot  = optField(opts, 'resultRoot',      '');
    configNames = optField(opts, 'configNames',     ...
        arrayfun(@(k) sprintf('cfg%d',k), (1:nConfigs)', 'UniformOutput', false));

    z       = optResult.zFinal(:);
    x_arm   = optResult.xFinal(:);
    C0      = optResult.C0(:);
    weights = weights(:);
    nElemsArm = analyses{1}.getTotalElemsNumber();
    assert(numel(x_arm) == nElemsArm, ...
        'postprocessSIMPResult: optResult.xFinal length %d does not match analysis element count %d.', ...
        numel(x_arm), nElemsArm);
    candidateReanalysis = runReanalysis || runSweep;

    %% -- Anchor elements (shared by all cleanup calls) -------------------
    anchorElems = findAnchorElements(analyses, model.mesh);

    %% -- 1. Volume-preserving threshold ----------------------------------
    fprintf('\n[postprocess] Volume-preserving threshold (targetVF=%.4f)...\n', VolFrac);
    t_vp     = findVolumeThreshold(x_arm, VolFrac);
    solid_vp_raw = x_arm >= t_vp;
    solid_vp = removeDisconnectedComponents(solid_vp_raw, model.mesh, anchorElems);

    vpCandidate = buildCandidate(solid_vp_raw, solid_vp, x_arm, analyses, penal, pAgg, weights, C0, ...
        candidateReanalysis, useParallel, sprintf('volume_t%.2f', t_vp));
    vpCandidate.threshold = t_vp;
    fprintf('  threshold=%.4f  V=%.4f  J_cont=%.6f\n', ...
        t_vp, vpCandidate.volFrac, vpCandidate.J_continuous);

    %% -- 2. Heaviside sharpening -----------------------------------------
    % Only applicable when the optimiser used a density filter + projection
    % (e.g. solveSIMPVolumeStressMMA with useProjection=true).  For
    % sensitivity-filtered SIMP (solveSIMPComplianceVolumeMMA) the design
    % variable IS the physical density, so applying Wfilter*z shifts
    % intermediate values and collapses volume — set opts.skipHeaviside=true
    % for those cases.
    heavisideCandidates = {};
    if ~skipHeaviside
        nBeta = numel(betaSeq);
        heavisideCandidates = cell(1, nBeta);

        for i = 1:nBeta
            beta = betaSeq(i);
            fprintf('[postprocess] Heaviside sharpening beta=%-4d...', beta);
            rho_sharp = heavisideSharpenPost(z, Wfilter, beta, eta, fixedVars);
            x_sharp   = designToArmDensity(rho_sharp, model, nElemsArm);

            solid_hi_raw = x_sharp >= 0.5;
            solid_hi  = removeDisconnectedComponents(solid_hi_raw, model.mesh, anchorElems);

            hiCandidate = buildCandidate(solid_hi_raw, solid_hi, x_sharp, analyses, penal, pAgg, weights, C0, ...
                candidateReanalysis, useParallel, sprintf('heaviside_b%d', beta));
            hiCandidate.beta = beta;
            heavisideCandidates{i} = hiCandidate;
            fprintf('  V=%.4f  J_cont=%.6f\n', hiCandidate.volFrac, hiCandidate.J_continuous);
        end
    end

    %% -- 3. Optional threshold sweep -------------------------------------
    sweepResult = [];
    sweepCandidate = [];
    if runSweep
        fprintf('[postprocess] Threshold sweep (%d points)...\n', numel(thresholds));
        sweepOpts = struct('pAgg', pAgg, 'weights', weights, 'C0', C0, ...
            'useParallel', useParallel, 'mesh', model.mesh, 'anchorElems', anchorElems);
        sweepResult = sweepThresholdsReanalysis(x_arm, analyses, thresholds, sweepOpts);
        t_sw = sweepResult.bestThreshold;
        solid_sw_raw = x_arm >= t_sw;
        solid_sw = removeDisconnectedComponents(solid_sw_raw, model.mesh, anchorElems);
        sweepCandidate = buildCandidate(solid_sw_raw, solid_sw, x_arm, analyses, penal, pAgg, weights, C0, ...
            true, useParallel, sprintf('sweep_t%.2f', t_sw));
        sweepCandidate.threshold = t_sw;
        sweepCandidate.J_binary = sweepResult.J(sweepResult.bestIdx);
        sweepCandidate.C_binary = sweepResult.C(sweepResult.bestIdx, :)';
        sweepCandidate.selectedScore = sweepCandidate.J_binary;
        sweepCandidate.selectedBy = "J_binary";
    end

    %% -- 4. Select best candidate ----------------------------------------
    allCandidates = [{vpCandidate}, heavisideCandidates];
    if ~isempty(sweepCandidate)
        allCandidates = [allCandidates, {sweepCandidate}];
    end
    [bestIdx, selectedBy] = selectBestComplianceCandidate(allCandidates, candidateReanalysis);
    bestCandidate = allCandidates{bestIdx};
    bestCandidate.selectedBy = selectedBy;
    fprintf('[postprocess] Best method: %s  (%s=%.6f, V=%.4f)\n', ...
        bestCandidate.label, char(selectedBy), bestCandidate.selectedScore, bestCandidate.volFrac);

    %% -- 5. Assemble output ----------------------------------------------
    postResult.volumeThreshold = vpCandidate;
    postResult.heaviside        = [heavisideCandidates{:}];
    postResult.sweepCandidate   = sweepCandidate;
    postResult.best             = bestCandidate;
    postResult.sweep            = sweepResult;

    %% -- 6. Save to disk -------------------------------------------------
    if ~isempty(resultRoot)
        savePostprocessResults(postResult, model, allCandidates, sweepResult, ...
            configNames, resultRoot);
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function cand = buildCandidate(solidRaw, solid, x_continuous, analyses, penal, pAgg, weights, C0, ...
        runReanalysis, useParallel, label)
% Evaluate a candidate topology.
% J_continuous: diagnostic SIMP objective using x_continuous.
% J_binary / C_binary: binary FE reanalysis (only if runReanalysis=true)

    nConfigs = numel(analyses);
    nElems = analyses{1}.getTotalElemsNumber();
    x_continuous = x_continuous(:);
    solidRaw = logical(solidRaw(:));
    solid = logical(solid(:));
    assert(numel(x_continuous) == nElems, ...
        'Candidate %s continuous density length %d does not match element count %d.', ...
        label, numel(x_continuous), nElems);
    assert(numel(solid) == nElems && numel(solidRaw) == nElems, ...
        'Candidate %s solid mask length does not match element count %d.', label, nElems);
    assert(all(isfinite(x_continuous)) && all(x_continuous >= 0 & x_continuous <= 1), ...
        'Candidate %s continuous density must be finite and in [0, 1].', label);

    xVoid = 1e-6;
    x_bin = double(solid) + (~solid) * xVoid;

    cand.label       = label;
    cand.xContinuous = x_continuous;
    cand.solidRaw    = solidRaw;
    cand.solid       = solid;
    cand.xBinary     = x_bin;
    cand.volFracRaw  = mean(solidRaw);
    cand.volFrac     = mean(solid);

    % Diagnostic objective on the continuous field.
    [J_s, ~, ~, C_s] = evaluateObjectiveAndGradient( ...
        analyses, x_continuous, penal, pAgg, weights, C0, useParallel);
    cand.J_continuous = J_s;
    cand.C_continuous = C_s;
    cand.J_simp = J_s; % Backward-compatible alias.
    cand.C_simp = C_s;

    % Structural performance metrics on the binary topology
    perfMetrics  = evaluateStructuralPerformance(analyses, x_bin, 1, useParallel);
    cand.sHM_max = perfMetrics.sHM_max;
    cand.u_max   = perfMetrics.u_max;

    % Binary FE reanalysis (expensive)
    cand.J_binary = NaN;
    cand.C_binary = NaN(nConfigs, 1);
    if runReanalysis
        C_bin  = zeros(nConfigs, 1);
        useP   = useParallel && license('test', 'Distrib_Computing_Toolbox');
        if useP
            parfor k = 1:nConfigs
                C_bin(k) = computeComplianceOnly(analyses{k}, x_bin, 1);
            end
        else
            for k = 1:nConfigs
                C_bin(k) = computeComplianceOnly(analyses{k}, x_bin, 1);
            end
        end
        Cagg           = C_bin ./ C0;
        weightedSum    = sum(weights .* (Cagg .^ pAgg));
        cand.J_binary  = weightedSum ^ (1.0 / pAgg);
        cand.C_binary  = C_bin;
    end
    if runReanalysis
        cand.selectedScore = cand.J_binary;
        cand.selectedBy = "J_binary";
    else
        cand.selectedScore = cand.J_continuous;
        cand.selectedBy = "J_continuous";
    end
end

% -------------------------------------------------------------------------
function [bestIdx, selectedBy] = selectBestComplianceCandidate(allCandidates, runReanalysis)
    if runReanalysis
        selectedBy = "J_binary";
        scores = cellfun(@(c) c.J_binary, allCandidates);
    else
        selectedBy = "J_continuous";
        scores = cellfun(@(c) c.J_continuous, allCandidates);
    end
    valid = isfinite(scores);
    assert(any(valid), 'postprocessSIMPResult: no finite %s values for candidate selection.', char(selectedBy));
    scores(~valid) = Inf;
    [~, bestIdx] = min(scores);
    best = allCandidates{bestIdx};
    best.selectedScore = scores(bestIdx);
    best.selectedBy = selectedBy;
    allCandidates{bestIdx} = best;
end

% -------------------------------------------------------------------------
function savePostprocessResults(~, model, allCandidates, sweepResult, ...
        ~, resultRoot)

    % --- Summary CSV ---
    nCands = numel(allCandidates);
    rows   = cell(nCands, 1);
    for i = 1:nCands
        c = allCandidates{i};
        row.label         = string(c.label);
        row.volFracRaw    = c.volFracRaw;
        row.volFracClean  = c.volFrac;
        row.J_continuous  = c.J_continuous;
        row.J_binary      = c.J_binary;
        row.selectedScore = c.selectedScore;
        row.selectedBy    = string(c.selectedBy);
        row.nSolidRaw     = sum(c.solidRaw);
        row.nSolidElems   = sum(c.solid);
        rows{i} = row;
    end
    T = struct2table(vertcat(rows{:}));
    writetable(T, fullfile(resultRoot, 'postprocess_summary.csv'));
    fprintf('[postprocess] Saved postprocess_summary.csv\n');

    % --- Topology plot per candidate ---
    for i = 1:nCands
        c   = allCandidates{i};
        stem = sprintf('postprocess_%s', c.label);
        fig  = figure('Visible', 'off', 'Name', c.label);
        hold on; axis off; daspect([1 1 1]); view(45, 35);
        model.fe.plotSolidSelected(model.mesh.nodes, c.solid, [0.45 0.60 0.80]);
        title(sprintf('%s  |  V=%.3f  %s=%.4f  sHM=%.3e  u=%.3e', ...
            strrep(c.label,'_',' '), c.volFrac, char(c.selectedBy), c.selectedScore, ...
            c.sHM_max, c.u_max), 'Interpreter', 'tex');
        saveas(fig, fullfile(resultRoot, [stem '.png']));
        set(fig, 'Visible', 'on');
        savefig(fig, fullfile(resultRoot, [stem '.fig']));
        close(fig);
    end

    % --- Pareto plot from sweep ---
    if ~isempty(sweepResult)
        fig = figure('Visible', 'off', 'Name', 'Threshold sweep Pareto');
        plot(sweepResult.volFrac, sweepResult.J, 'o-b', 'LineWidth', 1.5);
        hold on;
        plot(sweepResult.volFrac(sweepResult.bestIdx), sweepResult.J(sweepResult.bestIdx), ...
            'rs', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel('Volume fraction'); ylabel('Normalised J');
        title('Threshold sweep: J vs volume fraction');
        legend('sweep', 'best', 'Location', 'northwest');
        grid on;
        saveas(fig, fullfile(resultRoot, 'postprocess_pareto.png'));
        set(fig, 'Visible', 'on');
        savefig(fig, fullfile(resultRoot, 'postprocess_pareto.fig'));
        close(fig);
        writetable(struct2table(struct( ...
            'threshold', num2cell(sweepResult.thresholds), ...
            'volFrac',   num2cell(sweepResult.volFrac), ...
            'J',         num2cell(sweepResult.J))), ...
            fullfile(resultRoot, 'postprocess_sweep.csv'));
        fprintf('[postprocess] Saved postprocess_pareto.png and postprocess_sweep.csv\n');
    end
end

% -------------------------------------------------------------------------
function x_arm = designToArmDensity(x_design, model, nElemsArm)
    x_design = x_design(:);
    if numel(x_design) == nElemsArm
        x_arm = x_design;
    else
        x_arm = model.segmentToArm(x_design);
        assert(numel(x_arm) == nElemsArm, ...
            ['postprocessSIMPResult: mapped design length %d does not match ' ...
             'analysis element count %d.'], numel(x_arm), nElemsArm);
    end
end

% -------------------------------------------------------------------------
function v = optField(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

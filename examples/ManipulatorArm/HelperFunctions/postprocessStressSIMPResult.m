function postResult = postprocessStressSIMPResult(optResult, model, analyses, Wfilter, opts)
% POSTPROCESSSTRESSSIMPRESULT  Post-process a stress-constrained SIMP result.
%
%   postResult = postprocessStressSIMPResult(optResult, model, analyses, Wfilter, opts)
%
%   Applies density-threshold extraction strategies to the optimised density
%   field and selects the best topology by minimising volume fraction among
%   stress-feasible candidates.  When binary reanalysis is enabled,
%   feasibility is based on binary stress; otherwise continuous-field stress
%   is used as a diagnostic fallback.
%
%   Inputs:
%     optResult  struct from solveSIMPVolumeStressMMA.  Required fields:
%                  .zFinal          [H x 1]  raw (pre-filter) design variable,
%                                           or full-arm density for unlinked runs
%                  .xFinal          [nElems x 1]  full-arm physical density
%                  .finalVolumeFraction  scalar
%     model      ManipulatorModel3D
%     analyses   {nConfigs x 1}  FEAnalysis objects
%     Wfilter    [H x H]  sensitivity/density filter used during optimisation
%     opts       struct  (all fields optional unless marked required)
%
%   Required opts fields:
%     penal              SIMP exponent
%     stressPNorm        p-norm aggregation exponent
%     stressRelaxationQ  qp-relaxation exponent
%     nStressClusters    number of stress-level clusters per config
%     Starget            [nConfigs*nStressClusters x 1]  stress targets
%     constraintTol      feasibility threshold on g = S/Starget - 1
%
%   Optional opts fields:
%     betaSeq      Heaviside beta values to evaluate   (default: [8 16 32 64])
%     eta          Heaviside projection threshold       (default: 0.5)
%     fixedVars    indices forced to 1 in z-space       (default: [])
%     VolFrac      target for volume-preserving threshold
%                  (default: optResult.finalVolumeFraction)
%     runReanalysis  logical — re-evaluate stress on binary field
%                    (expensive; default: false)
%     useParallel  logical                              (default: false)
%     resultRoot   path for saving outputs              (default: '' = no saving)
%     configNames  {nConfigs x 1} string labels
%
%   Output: postResult struct with fields
%     volumeThreshold   struct  volume-preserving extraction candidate
%     heaviside(i)      struct  Heaviside at betaSeq(i)
%     best              struct  copy of selected best candidate
%
%   Each candidate struct has:
%     .xContinuous        [nElems x 1] double   continuous density used diagnostically
%     .solidRaw           [nElems x 1] logical  threshold mask before cleanup
%     .solid              [nElems x 1] logical  binary mask after cleanup
%     .xBinary            [nElems x 1] double   solid + xVoid*(~solid)
%     .volFracRaw         scalar
%     .volFrac            scalar
%     .S_continuous       [nConstraints x 1]  p-norm stress at sharpened field
%     .maxStress_continuous [nConstraints x 1]
%     .maxConstraintContinuous  scalar  max(S/Starget - 1)
%     .feasibleContinuous logical  (maxConstraintContinuous <= constraintTol)
%     .S_binary           [nConstraints x 1]  (if runReanalysis=true, else NaN)
%     .maxConstraintBinary  scalar  (if runReanalysis=true, else NaN)
%     .feasibleSelected   logical
%     .selectedScore      scalar
%     .selectedBy         string
%     .label              string

    %% -- Validate required opts ---------------------------------------------
    required = {'penal','stressPNorm','stressRelaxationQ','nStressClusters','Starget','constraintTol'};
    for i = 1:numel(required)
        assert(isfield(opts, required{i}), ...
            'postprocessStressSIMPResult: opts.%s is required.', required{i});
    end

    nConfigs        = numel(analyses);
    penal           = opts.penal;
    stressPNorm     = opts.stressPNorm;
    q               = opts.stressRelaxationQ;
    nStressClusters = opts.nStressClusters;
    Starget         = opts.Starget(:);
    constraintTol   = opts.constraintTol;

    betaSeq       = optField(opts, 'betaSeq',       [8 16 32 64]);
    eta           = optField(opts, 'eta',           0.5);
    fixedVars     = optField(opts, 'fixedVars',     []);
    VolFrac       = optField(opts, 'VolFrac',       optResult.finalVolumeFraction);
    runReanalysis = optField(opts, 'runReanalysis', false);
    useParallel   = optField(opts, 'useParallel',   false);
    resultRoot    = optField(opts, 'resultRoot',    '');
    configNames   = optField(opts, 'configNames',   ...
        arrayfun(@(k) sprintf('cfg%d',k), (1:nConfigs)', 'UniformOutput', false));

    z     = optResult.zFinal(:);
    x_arm = optResult.xFinal(:);
    nElemsArm = analyses{1}.getTotalElemsNumber();
    assert(numel(x_arm) == nElemsArm, ...
        'postprocessStressSIMPResult: optResult.xFinal length %d does not match analysis element count %d.', ...
        numel(x_arm), nElemsArm);

    %% -- Anchor elements ----------------------------------------------------
    anchorElems = findAnchorElements(analyses, model.mesh);

    %% -- 1. Volume-preserving threshold ------------------------------------
    fprintf('\n[postprocess stress] Volume-preserving threshold (targetVF=%.4f)...\n', VolFrac);
    t_vp     = findVolumeThreshold(x_arm, VolFrac);
    solid_vp_raw = x_arm >= t_vp;
    solid_vp = removeDisconnectedComponents(solid_vp_raw, model.mesh, anchorElems);

    vpCandidate = buildCandidateStress(solid_vp_raw, solid_vp, x_arm, analyses, penal, stressPNorm, q, ...
        Starget, nStressClusters, constraintTol, runReanalysis, useParallel, ...
        sprintf('volume_t%.2f', t_vp));
    vpCandidate.threshold = t_vp;
    fprintf('  threshold=%.4f  V=%.4f  maxG_cont=%.4f  feasible=%d\n', ...
        t_vp, vpCandidate.volFrac, vpCandidate.maxConstraintContinuous, vpCandidate.feasibleSelected);

    %% -- 2. Heaviside sharpening -------------------------------------------
    nBeta = numel(betaSeq);
    heavisideCandidates = cell(1, nBeta);

    for i = 1:nBeta
        beta = betaSeq(i);
        fprintf('[postprocess stress] Heaviside sharpening beta=%-4d...', beta);
        rho_sharp = heavisideSharpenPost(z, Wfilter, beta, eta, fixedVars);
        x_sharp   = designToArmDensity(rho_sharp, model, nElemsArm);

        solid_hi_raw = x_sharp >= 0.5;
        solid_hi  = removeDisconnectedComponents(solid_hi_raw, model.mesh, anchorElems);

        hiCandidate = buildCandidateStress(solid_hi_raw, solid_hi, x_sharp, analyses, penal, stressPNorm, q, ...
            Starget, nStressClusters, constraintTol, runReanalysis, useParallel, ...
            sprintf('heaviside_b%d', beta));
        hiCandidate.beta = beta;
        heavisideCandidates{i} = hiCandidate;
        fprintf('  V=%.4f  maxG_cont=%.4f  feasible=%d\n', ...
            hiCandidate.volFrac, hiCandidate.maxConstraintContinuous, hiCandidate.feasibleSelected);
    end

    %% -- 3. Select best candidate ------------------------------------------
    allCandidates = [{vpCandidate}, heavisideCandidates];
    [bestIdx, selectedBy] = selectBestStressCandidate(allCandidates, runReanalysis);

    bestCandidate = allCandidates{bestIdx};
    bestCandidate.selectedBy = selectedBy;
    fprintf('[postprocess stress] Best: %s  (V=%.4f, %s=%.4f, feasible=%d)\n', ...
        bestCandidate.label, bestCandidate.volFrac, ...
        char(selectedBy), bestCandidate.selectedScore, bestCandidate.feasibleSelected);

    %% -- 4. Assemble output ------------------------------------------------
    postResult.volumeThreshold = vpCandidate;
    postResult.heaviside        = [heavisideCandidates{:}];
    postResult.best             = bestCandidate;

    %% -- 5. Save to disk ---------------------------------------------------
    if ~isempty(resultRoot)
        saveStressPostprocessResults(model, allCandidates, configNames, resultRoot);
    end
end

% =========================================================================
% Local helpers
% =========================================================================

function cand = buildCandidateStress(solidRaw, solid, x_continuous, analyses, penal, stressPNorm, q, ...
        Starget, nStressClusters, constraintTol, runReanalysis, useParallel, label)

    nConstraints = numel(Starget);
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

    % Fast stress evaluation on the continuous sharpened field
    [S_cont, ms_cont] = evaluateStressSetOnly(analyses, x_continuous, penal, stressPNorm, ...
        q, Starget, useParallel, nStressClusters);
    cand.S_continuous           = S_cont;
    cand.maxStress_continuous   = ms_cont;
    cand.maxConstraintContinuous = max(S_cont ./ max(Starget(:), eps) - 1.0);
    cand.feasibleContinuous      = cand.maxConstraintContinuous <= constraintTol;
    cand.feasible                = cand.feasibleContinuous; % Backward-compatible alias.

    % Structural performance metrics on the binary topology
    perfMetrics  = evaluateStructuralPerformance(analyses, x_bin, 1, useParallel);
    cand.sHM_max = perfMetrics.sHM_max;
    cand.u_max   = perfMetrics.u_max;

    % Binary reanalysis (expensive, opt-in)
    cand.S_binary            = NaN(nConstraints, 1);
    cand.maxConstraintBinary = NaN;
    if runReanalysis
        [S_bin, ~] = evaluateStressSetOnly(analyses, x_bin, 1, stressPNorm, q, ...
            Starget, useParallel, nStressClusters);
        cand.S_binary            = S_bin;
        cand.maxConstraintBinary = max(S_bin ./ max(Starget(:), eps) - 1.0);
    end
    if runReanalysis
        cand.selectedScore = cand.maxConstraintBinary;
        cand.selectedBy = "maxConstraintBinary";
        cand.feasibleSelected = cand.maxConstraintBinary <= constraintTol;
    else
        cand.selectedScore = cand.maxConstraintContinuous;
        cand.selectedBy = "maxConstraintContinuous";
        cand.feasibleSelected = cand.feasibleContinuous;
    end
end

% -------------------------------------------------------------------------
function [bestIdx, selectedBy] = selectBestStressCandidate(allCandidates, runReanalysis)
    if runReanalysis
        selectedBy = "maxConstraintBinary";
        scores = cellfun(@(c) c.maxConstraintBinary, allCandidates);
        feasible = cellfun(@(c) c.feasibleSelected, allCandidates);
    else
        selectedBy = "maxConstraintContinuous";
        scores = cellfun(@(c) c.maxConstraintContinuous, allCandidates);
        feasible = cellfun(@(c) c.feasibleSelected, allCandidates);
    end
    valid = isfinite(scores);
    assert(any(valid), 'postprocessStressSIMPResult: no finite %s values for candidate selection.', char(selectedBy));

    if any(feasible & valid)
        vols = cellfun(@(c) c.volFrac, allCandidates);
        vols(~(feasible & valid)) = Inf;
        [~, bestIdx] = min(vols);
    else
        scores(~valid) = Inf;
        [~, bestIdx] = min(scores);
        fprintf('[postprocess stress] Warning: no feasible candidate found by %s.\n', char(selectedBy));
    end
end

% -------------------------------------------------------------------------
function saveStressPostprocessResults(model, allCandidates, configNames, resultRoot)

    % Summary CSV
    nCands = numel(allCandidates);
    rows   = cell(nCands, 1);
    for i = 1:nCands
        c = allCandidates{i};
        row.label                = string(c.label);
        row.volFracRaw           = c.volFracRaw;
        row.volFracClean         = c.volFrac;
        row.maxConstraintCont    = c.maxConstraintContinuous;
        row.feasibleContinuous   = double(c.feasibleContinuous);
        row.maxConstraintBinary  = c.maxConstraintBinary;
        row.feasibleSelected     = double(c.feasibleSelected);
        row.selectedScore        = c.selectedScore;
        row.selectedBy           = string(c.selectedBy);
        row.nSolidRaw            = sum(c.solidRaw);
        row.nSolidElems          = sum(c.solid);
        rows{i} = row;
    end
    T = struct2table(vertcat(rows{:}));
    writetable(T, fullfile(resultRoot, 'postprocess_stress_summary.csv'));
    fprintf('[postprocess stress] Saved postprocess_stress_summary.csv\n');

    % Topology plot per candidate
    for i = 1:nCands
        c    = allCandidates{i};
        stem = sprintf('postprocess_stress_%s', c.label);
        fig  = figure('Visible', 'off', 'Name', c.label);
        hold on; axis off; daspect([1 1 1]); view(45, 35);
        model.fe.plotSolidSelected(model.mesh.nodes, c.solid, [0.45 0.60 0.80]);
        title(sprintf('%s  |  V=%.3f  %s=%.3f  sHM=%.3e  u=%.3e', ...
            strrep(c.label,'_',' '), c.volFrac, char(c.selectedBy), c.selectedScore, ...
            c.sHM_max, c.u_max), 'Interpreter', 'tex');
        saveas(fig, fullfile(resultRoot, [stem '.png']));
        set(fig, 'Visible', 'on');
        savefig(fig, fullfile(resultRoot, [stem '.fig']));
        close(fig);
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
            ['postprocessStressSIMPResult: mapped design length %d does not match ' ...
             'analysis element count %d.'], numel(x_arm), nElemsArm);
    end
end

% -------------------------------------------------------------------------
function v = optField(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

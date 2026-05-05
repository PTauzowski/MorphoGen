function result = evaluateStructuralPerformance(analyses, x, penal, useParallel)
% EVALUATESTRUCTURALPERFORMANCE  Max HM stress and max displacement for a density field.
%
%   result = evaluateStructuralPerformance(analyses, x, penal, useParallel)
%
%   For each load configuration, solves the weighted FEM problem and extracts:
%     - Max element-averaged nodal Huber-Mises (HM) stress
%     - Max nodal displacement magnitude (L2 norm of ux,uy,uz)
%
%   Inputs:
%     analyses     {nConfigs x 1}  FEAnalysis objects
%     x            [nElems x 1]   physical density field
%                  Pass ones(nElems,1) with penal=1 for the full-material reference.
%     penal        SIMP penalisation exponent (1 = linear, 3 = standard SIMP)
%     useParallel  logical  (default false)
%
%   Output: result struct with fields
%     .sHM_perConfig  [nConfigs x 1]  max HM stress per config
%     .u_perConfig    [nConfigs x 1]  max displacement magnitude per config
%     .sHM_max        scalar          max over all configs
%     .u_max          scalar          max over all configs

    if nargin < 4, useParallel = false; end

    nConfigs    = numel(analyses);
    x           = max(x(:), 1e-9);
    xPenal      = x .^ penal;
    useParallel = useParallel && license('test', 'Distrib_Computing_Toolbox');

    sHM_k = zeros(nConfigs, 1);
    u_k   = zeros(nConfigs, 1);

    if useParallel
        parfor k = 1:nConfigs
            [sHM_k(k), u_k(k)] = evalOneConfig(analyses{k}, xPenal);
        end
    else
        for k = 1:nConfigs
            [sHM_k(k), u_k(k)] = evalOneConfig(analyses{k}, xPenal);
        end
    end

    result.sHM_perConfig = sHM_k;
    result.u_perConfig   = u_k;
    result.sHM_max       = max(sHM_k);
    result.u_max         = max(u_k);
end

% =========================================================================
% Local helpers
% =========================================================================

function [sHM_max, u_max] = evalOneConfig(analysis, xPenal)
    nElems = analysis.getTotalElemsNumber();

    analysis.solveWeighted(xPenal, false);
    analysis.computeElementResults(xPenal);

    % Max Huber-Mises stress: mean of nodal sHM over each element's nodes
    elemIndices = analysis.getElemIndices();
    sHM = zeros(nElems, 1);
    for i = 1:numel(analysis.felems)
        fe     = analysis.felems{i};
        hmIdx  = find(fe.results.names == "sHM", 1);
        if isempty(hmIdx), continue; end
        ids    = elemIndices{i};
        for j = 1:numel(ids)
            nodeIds     = fe.elems(j, :);
            sHM(ids(j)) = mean(fe.results.nodal.all(nodeIds, hmIdx));
        end
    end
    sHM_max = max(sHM);

    % Max displacement magnitude (L2 norm of ux,uy,uz at each node)
    dispIdx = analysis.findDOFsIndices(["ux", "uy", "uz"]);
    u_max   = max(vecnorm(analysis.qnodal(:, dispIdx), 2, 2));
end

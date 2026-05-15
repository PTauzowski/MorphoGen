function writeTopologyPerformanceComparisonCsv(resultRoot, fileName, fullPipeMetrics, candidates, bestLabel)
% writeTopologyPerformanceComparisonCsv  Compare extracted topologies to full pipe.

    if nargin < 5
        bestLabel = "";
    end
    if ~iscell(candidates)
        candidates = num2cell(candidates);
    end

    n = numel(candidates);
    rows = cell(n, 1);
    for i = 1:n
        c = candidates{i};
        row.label = string(c.label);
        row.is_best = string(c.label) == string(bestLabel);
        row.vol_frac = valueOrNaN(c, 'volFrac');
        row.n_solid_elems = valueOrNaN(c, 'nSolidElems');

        row.sHM_full_pipe = fullPipeMetrics.sHM_max;
        row.sHM_topology = valueOrNaN(c, 'sHM_max');
        row.sHM_ratio_topology_to_full_pipe = row.sHM_topology / max(row.sHM_full_pipe, eps);

        row.u_full_pipe = fullPipeMetrics.u_max;
        row.u_topology = valueOrNaN(c, 'u_max');
        row.u_ratio_topology_to_full_pipe = row.u_topology / max(row.u_full_pipe, eps);

        row.uz_full_pipe = fullPipeMetrics.uz_max;
        row.uz_topology = valueOrNaN(c, 'uz_max');
        row.uz_ratio_topology_to_full_pipe = row.uz_topology / max(row.uz_full_pipe, eps);
        rows{i} = row;
    end

    T = struct2table(vertcat(rows{:}));
    writetable(T, fullfile(resultRoot, fileName));
    if n > 0
        bestRows = T(T.is_best, :);
        if isempty(bestRows)
            bestRows = T(1, :);
        end
        fprintf('[topology comparison] full pipe -> %s: sHM %.3e -> %.3e (x%.2f), uz %.3e -> %.3e (x%.2f)\n', ...
            char(bestRows.label(1)), fullPipeMetrics.sHM_max, bestRows.sHM_topology(1), ...
            bestRows.sHM_ratio_topology_to_full_pipe(1), fullPipeMetrics.uz_max, ...
            bestRows.uz_topology(1), bestRows.uz_ratio_topology_to_full_pipe(1));
    end
end

function v = valueOrNaN(s, field)
    if isfield(s, field)
        v = s.(field);
    elseif strcmp(field, 'nSolidElems') && isfield(s, 'solid')
        v = sum(s.solid);
    else
        v = NaN;
    end
end

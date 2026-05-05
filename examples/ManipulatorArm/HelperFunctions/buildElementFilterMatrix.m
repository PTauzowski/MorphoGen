function W = buildElementFilterMatrix(nodes, elems, elemIds, Rfilter)
% buildElementFilterMatrix  Row-normalized element-centroid density filter.
%
%   W = buildElementFilterMatrix(nodes, elems, elemIds, Rfilter) builds a
%   sparse matrix so filteredSensitivity = W * sensitivity.  elemIds selects
%   the design elements to include in the filter space.

    elemIds = elemIds(:);
    n = numel(elemIds);
    centroids = zeros(n, size(nodes, 2));
    for i = 1:n
        centroids(i, :) = mean(nodes(elems(elemIds(i), :), :), 1);
    end

    ii = cell(n, 1);
    jj = cell(n, 1);
    vv = cell(n, 1);

    for i = 1:n
        distances = vecnorm(centroids - centroids(i, :), 2, 2);
        neighbours = find(distances <= Rfilter);
        weights = Rfilter - distances(neighbours);
        weights = weights / max(sum(weights), eps);

        ii{i} = repmat(i, numel(neighbours), 1);
        jj{i} = neighbours(:);
        vv{i} = weights(:);
    end

    W = sparse(vertcat(ii{:}), vertcat(jj{:}), vertcat(vv{:}), n, n);
end

function rhoBinary = makeVolumePreservingBinaryDensity(rho, targetVf, fixedIds)
% makeVolumePreservingBinaryDensity  Keep the densest free elements at target volume.

    rho = rho(:);
    n = numel(rho);
    fixed = false(n, 1);
    fixed(fixedIds(:)) = true;

    rhoBinary = zeros(n, 1);
    rhoBinary(fixed) = 1.0;

    nTarget = round(targetVf * n);
    nFreeKeep = max(0, min(nnz(~fixed), nTarget - nnz(fixed)));
    freeIds = find(~fixed);
    if nFreeKeep > 0
        [~, order] = sort(rho(freeIds), 'descend');
        rhoBinary(freeIds(order(1:nFreeKeep))) = 1.0;
    end
end

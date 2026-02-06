function [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
    edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, eigTol)
%LOWRANKFROMELEMENTUPDATES_VECKE_ARRAY Build Woodbury factors from element updates.
%
%   [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
%       edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, eigTol)
%
% Inputs:
%   edofs_list   : (nChanged x ndLocal) global DOFs for changed elements.
%   Ke0_vec_list : either cell(nChanged,1) of (ndLocal^2)x1 vectors or
%                  numeric matrix (nChanged x ndLocal^2) / (ndLocal^2 x nChanged).
%   w_old,w_new  : (nChanged x 1) old/new element stiffness weights.
%   freeMap      : (nGlobalDof x 1) map global DOF -> free index (0 = fixed).
%   nFree        : number of free DOFs.
%   eigTol       : relative eigenvalue truncation tolerance.
%
% Outputs:
%   U            : (nFree x m) low-rank basis.
%   C            : (m x m) diagonal sign matrix (+1/-1 per retained mode).
%   info         : diagnostics struct with fields:
%                  rank, nChanged, keptPerElem, eigTol, nZeroDw, nNoFree,
%                  maxSymErr, nNonSym.
%
% Notes:
%   - Does not assemble a global dK.
%   - For each changed element i:
%       dKe_i = (w_new(i)-w_old(i)) * Ke0_i
%     is restricted to free local DOFs and eigen-truncated.
%   - The decomposition dK ~= U*C*U' is exact up to eigTol truncation.
%   - Works for H8 (ndLocal=24), Q4 (ndLocal=8), and other fixed-size elements.

    if nargin < 7 || isempty(eigTol)
        eigTol = 1.0e-10;
    end

    assert(isnumeric(edofs_list) && ismatrix(edofs_list), ...
        'edofs_list must be a numeric matrix.');
    ndLocal = size(edofs_list, 2);
    assert(ndLocal > 0, 'edofs_list must have at least one local DOF column.');
    expectedKeLen = ndLocal * ndLocal;
    assert(isscalar(nFree) && nFree >= 0 && nFree == floor(nFree), ...
        'nFree must be a non-negative integer scalar.');
    assert(isvector(freeMap), 'freeMap must be a vector.');
    assert(isscalar(eigTol) && eigTol >= 0, ...
        'eigTol must be a non-negative scalar.');

    nChanged = size(edofs_list, 1);
    w_old = w_old(:);
    w_new = w_new(:);
    freeMap = freeMap(:);

    assert(numel(w_old) == nChanged, ...
        'w_old length must match size(edofs_list,1).');
    assert(numel(w_new) == nChanged, ...
        'w_new length must match size(edofs_list,1).');

    if nChanged > 0
        assert(all(edofs_list(:) >= 1) && all(edofs_list(:) == floor(edofs_list(:))), ...
            'edofs_list must contain positive integer DOF IDs.');
        assert(max(edofs_list(:)) <= numel(freeMap), ...
            'freeMap is too short for edofs_list entries.');
    end

    [fetchKeVec, keStorageTag] = localKeAccessor(Ke0_vec_list, nChanged, expectedKeLen);

    keptPerElem = zeros(nChanged, 1);
    nZeroDw = 0;
    nNoFree = 0;
    nNonSym = 0;
    maxSymErr = 0.0;
    m = 0;

    for i = 1:nChanged
        dw = w_new(i) - w_old(i);
        if dw == 0
            nZeroDw = nZeroDw + 1;
            continue;
        end

        Kevec = fetchKeVec(i);
        assert(numel(Kevec) == expectedKeLen, ...
            'Each Ke0 vector must have length ndLocal^2.');

        dKe = reshape(dw * Kevec(:), ndLocal, ndLocal);
        symErr = norm(dKe - dKe.', 'fro');
        maxSymErr = max(maxSymErr, symErr);
        if symErr > 1.0e-10 * max(1.0, norm(dKe, 'fro'))
            nNonSym = nNonSym + 1;
        end
        dKe = 0.5 * (dKe + dKe.');

        edofs = edofs_list(i, :).';
        fidx = freeMap(edofs);
        keep = fidx > 0;

        if ~any(keep)
            nNoFree = nNoFree + 1;
            continue;
        end

        fidx = fidx(keep);
        assert(all(fidx <= nFree), 'freeMap contains index above nFree.');
        loc = find(keep); %#ok<FNDSB>
        dKeF = dKe(loc, loc);
        dKeF = 0.5 * (dKeF + dKeF.');

        [~, D] = eig(dKeF);
        lam = real(diag(D));
        lamMax = max(abs(lam));
        thr = eigTol * max(1.0, lamMax);
        nKeep = nnz(abs(lam) > thr);

        keptPerElem(i) = nKeep;
        m = m + nKeep;
    end

    if nNonSym > 0
        warning('lowRankFromElementUpdates_vecKe_array:NonSymmetricUpdate', ...
            ['Detected %d non-symmetric local updates (max ||dKe-dKe''||_F = %.3e). ', ...
             'Updates were symmetrized. Input Ke storage: %s.'], ...
            nNonSym, maxSymErr, keStorageTag);
    end

    if m == 0
        U = zeros(nFree, 0);
        C = zeros(0, 0);
        info = struct( ...
            'rank', 0, ...
            'nChanged', nChanged, ...
            'ndLocal', ndLocal, ...
            'keptPerElem', keptPerElem, ...
            'eigTol', eigTol, ...
            'nZeroDw', nZeroDw, ...
            'nNoFree', nNoFree, ...
            'maxSymErr', maxSymErr, ...
            'nNonSym', nNonSym);
        return;
    end

    U = zeros(nFree, m);
    Cvals = zeros(m, 1);
    col = 1;

    for i = 1:nChanged
        dw = w_new(i) - w_old(i);
        if dw == 0
            continue;
        end

        Kevec = fetchKeVec(i);
        dKe = reshape(dw * Kevec(:), ndLocal, ndLocal);
        dKe = 0.5 * (dKe + dKe.');

        edofs = edofs_list(i, :).';
        fidx = freeMap(edofs);
        keep = fidx > 0;
        if ~any(keep)
            continue;
        end

        fidx = fidx(keep);
        loc = find(keep); %#ok<FNDSB>
        dKeF = dKe(loc, loc);
        dKeF = 0.5 * (dKeF + dKeF.');

        [V, D] = eig(dKeF);
        lam = real(diag(D));
        V = real(V);

        lamMax = max(abs(lam));
        thr = eigTol * max(1.0, lamMax);
        sel = find(abs(lam) > thr); %#ok<FNDSB>
        if isempty(sel)
            continue;
        end

        for j = 1:numel(sel)
            lambda = lam(sel(j));
            uj = zeros(nFree, 1);
            uj(fidx) = sqrt(abs(lambda)) * V(:, sel(j));

            U(:, col) = uj;
            Cvals(col) = sign(lambda);
            col = col + 1;
        end
    end

    if col <= m
        U = U(:, 1:col-1);
        Cvals = Cvals(1:col-1);
    end

    C = diag(Cvals);
    info = struct( ...
        'rank', size(U, 2), ...
        'nChanged', nChanged, ...
        'ndLocal', ndLocal, ...
        'keptPerElem', keptPerElem, ...
        'eigTol', eigTol, ...
        'nZeroDw', nZeroDw, ...
        'nNoFree', nNoFree, ...
        'maxSymErr', maxSymErr, ...
        'nNonSym', nNonSym);
end

function [fetchKeVec, storageTag] = localKeAccessor(Ke0_vec_list, nChanged, expectedKeLen)
    if iscell(Ke0_vec_list)
        assert(numel(Ke0_vec_list) == nChanged, ...
            'Ke0_vec_list cell length must match nChanged.');
        fetchKeVec = @(i) Ke0_vec_list{i};
        storageTag = 'cell';
        return;
    end

    assert(isnumeric(Ke0_vec_list) && ismatrix(Ke0_vec_list), ...
        'Ke0_vec_list must be a cell array or numeric matrix.');

    if size(Ke0_vec_list, 1) == nChanged && size(Ke0_vec_list, 2) == expectedKeLen
        fetchKeVec = @(i) Ke0_vec_list(i, :).';
        storageTag = sprintf('matrix(nChanged x %d)', expectedKeLen);
        return;
    end

    if size(Ke0_vec_list, 1) == expectedKeLen && size(Ke0_vec_list, 2) == nChanged
        fetchKeVec = @(i) Ke0_vec_list(:, i);
        storageTag = sprintf('matrix(%d x nChanged)', expectedKeLen);
        return;
    end

    error(['Ke0_vec_list must be cell(nChanged,1), (nChanged x ndLocal^2), ', ...
           'or (ndLocal^2 x nChanged).']);
end

% WoodburyReanalysisDemoH8
% Unit-like demo for H8 low-rank reanalysis without assembling global dK.

clear;
clc;
rng(2);

% --- Synthetic global system (free + constrained DOFs)
nGlobal = 180;
supports = false(nGlobal, 1);
supports(1:30) = true;

[II, JJ] = ndgrid(1:nGlobal, 1:nGlobal);
I = II(:);
J = JJ(:);
solver = LinearEquationsSystem(I, J, supports);

freeMap = solver.createFreeMap();
nFree = numel(solver.freedofs);

% Strong SPD base matrix on free DOFs.
A = randn(nFree, nFree);
Kbase_ff = A.' * A + 200.0 * speye(nFree);
Kbase = zeros(nGlobal, nGlobal);
Kbase(solver.freedofs, solver.freedofs) = Kbase_ff;
Kbase = 0.5 * (Kbase + Kbase.');
Kbase_vec = Kbase(:);
decBase = solver.setBaseMatrix(Kbase_vec);

% --- Synthetic "changed H8 elements" with local 24x24 updates
nChanged = 5;
edofs_list = zeros(nChanged, 24);
Ke0_vec_list = zeros(nChanged, 576);
w_old = 0.4 + 0.4 * rand(nChanged, 1);
w_new = w_old + 0.08 * randn(nChanged, 1); % allows both + and - updates

for e = 1:nChanged
    edofs_list(e, :) = randperm(nGlobal, 24);
    B = randn(24, 24);
    Ke0 = B.' * B + 1.0e-3 * eye(24);
    Ke0_vec_list(e, :) = Ke0(:).';
end

% Build reference current matrix K = Kbase + sum(dKe_e) for verification only.
Kcurr = Kbase;
for e = 1:nChanged
    dw = w_new(e) - w_old(e);
    dKe = reshape(dw * Ke0_vec_list(e, :).', 24, 24);
    dKe = 0.5 * (dKe + dKe.');
    idx = edofs_list(e, :);
    Kcurr(idx, idx) = Kcurr(idx, idx) + dKe;
end
Kcurr = 0.5 * (Kcurr + Kcurr.');
Kcurr_vec = Kcurr(:);

P = zeros(nGlobal, 1);
P(solver.freedofs) = randn(nFree, 1);

eigTol = 1.0e-10;
[U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
    edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, eigTol);

mMax = 120;
[qWoodbury, wbLog] = solver.solveWoodburyWithFallback( ...
    decBase, Kcurr_vec, P, U, C, info, mMax, 'direct');
qDirect = solver.solve(Kcurr_vec, P);

relErr = norm(qWoodbury - qDirect) / max(1.0, norm(qDirect));
fprintf('Woodbury rank=%d, method=%s, relative error=%.3e\n', ...
    info.rank, wbLog.method, relErr);
assert(relErr < 1.0e-8, ...
    'Woodbury reanalysis mismatch: relative error too large.');

% Force fallback branch for a quick consistency check.
[qFallback, fbLog] = solver.solveWoodburyWithFallback( ...
    decBase, Kcurr_vec, P, U, C, info, 1, 'direct');
relErrFallback = norm(qFallback - qDirect) / max(1.0, norm(qDirect));
fprintf('Fallback method=%s, relative error=%.3e\n', ...
    fbLog.method, relErrFallback);
assert(relErrFallback < 1.0e-12, ...
    'Fallback solve mismatch: relative error too large.');

% -------------------------------------------------------------------------
% TO-loop call pattern (H8, ndof/node = 3):
%
% changed = find(abs(w_new_all - w_old_all) > 0);
% if ~isempty(changed)
%     allEdofs = reshape( ...
%         (repelem(felem.elems, 1, 3) - 1) * 3 + repmat(1:3, size(felem.elems, 1), 8), ...
%         size(felem.elems, 1), 24);
%     edofs_list = allEdofs(changed, :);
%     w_old = w_old_all(changed);
%     w_new = w_new_all(changed);
%     Ke0_vec_list = Ke0_bank(changed, :);  % (nChanged x 576) or cell
%
%     [U, C, info] = lowRankFromElementUpdates_vecKe_array( ...
%         edofs_list, Ke0_vec_list, w_old, w_new, freeMap, nFree, eigTol);
%
%     [qfem, wbLog] = solver.solveWoodburyWithFallback( ...
%         decBase, Kcurr_vec, Pfem, U, C, info, mMax, 'direct');
% else
%     qfem = solver.solveWoodbury(decBase, Pfem, [], []);
% end
% -------------------------------------------------------------------------

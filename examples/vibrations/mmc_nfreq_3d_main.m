function results = mmc_nfreq_3d_main(config)
% MMC 3D natural frequency optimization (refactored driver)
% - Mirrors Appendix G logic with cleaned utilities.
% - Leaves 2D solver untouched.
%
% Call selftest_3d to run a quick 3-iteration sanity check.

if nargin < 1, config = struct(); end
cfg = defaults3d();
cfg = override(cfg, config);

[KE, ME] = fe_hex8(cfg.E, cfg.nu, cfg.rho, cfg.EL, cfg.EW, cfg.EH);
grid = grid_structured_3d(cfg.nelx, cfg.nely, cfg.nelz, cfg.EL, cfg.EW, cfg.EH);
bc   = bc_builders_3d(cfg.bcCase, cfg.nelx, cfg.nely, cfg.nelz);

nDsvb = length(cfg.dd(:));
nEhcp = nDsvb / cfg.N;

fsparse = @(i,j,s,sz) sparse(i,j,s,sz(1),sz(2));

xval  = cfg.dd(:); xold1 = xval; xold2 = xval;
xmin  = repmat(cfg.xmin, cfg.N, 1); xmax = repmat(cfg.xmax, cfg.N, 1);
low   = xmin; upp = xmax;

iter = 1; OBJ = []; CONS = []; objVr5 = 1; loop = []; var_obj = [];
actComp = 1:cfg.N; actDsvb = 1:nDsvb;
allPhi = [zeros(cfg.nNod, cfg.N) cfg.PhiNd];
allPhiDrv = sparse(cfg.nNod, nDsvb);

while objVr5 > 1e-4 && iter <= cfg.maxiter
    Time_iteration = clock;
    epsilon = continuation_eps(iter, cfg.ep_old, cfg.rank, cfg.mid_num, cfg.ini);

    % TDFs
    for i = actComp
        [allPhi, allPhiDrv, xval, actComp, actDsvb] = calc_Phi_3d(allPhi, allPhiDrv, xval, i, grid.LSgrid, cfg.p, nEhcp, epsilon, actComp, actDsvb, cfg.minSz);
    end
    allPhiAct = [allPhi(:, actComp) cfg.PhiNd];
    temp = exp(cfg.lmd * allPhiAct);
    Phimax = max(-1e3, log(sum(temp, 2)) / cfg.lmd);
    allPhiDrvAct = allPhiDrv(:, actDsvb);
    Phimaxdphi = kron(temp(:, 1:length(actComp)) ./ (sum(temp, 2) + eps), ones(1, nEhcp));
    PhimaxDrvAct = Phimaxdphi .* allPhiDrvAct;

    % densities
    H = heaviside_smooth(Phimax, cfg.alpha, epsilon);
    den = sum(H(grid.eleNodesID), 2) / 8;

    % connectivity check
    bw_den = reshape(den, cfg.nelx, cfg.nely, cfg.nelz) > cfg.alpha;
    bw_den_sol = reshape(bwlabeln(bw_den), cfg.nelx*cfg.nely*cfg.nelz, 1);
    bw_fix  = unique(bw_den_sol(bc.fixEle));
    bw_mass = nonzeros(unique(bw_den_sol(bc.massEle)));

    denSld = den; edofMatLft = grid.edofMat; freedofLft = bc.freeDof; Iar = grid.Iar0;
    if sum(ismember(bw_fix, bw_mass)) > 0
        eleLft = find(bw_den_sol == bw_mass);
        denSld = zeros(size(den)); denSld(eleLft) = den(eleLft);
        edofMatLft = grid.edofMat(eleLft, :);
        freedofLft = setdiff(edofMatLft(:), bc.fixDof);
        [iK1, jK1] = deal(edofMatLft(:, grid.sI)', edofMatLft(:, grid.sII)');
        Iar = sort([iK1(:), jK1(:)], 2, 'descend');
    end

    % assemble
    sK = reshape(KE(:) * denSld(:)', length(KE)*length(denSld), 1);
    K = fsparse(Iar(:,1), Iar(:,2), sK, [cfg.nDof, cfg.nDof]);
    K = K + K' - diag(diag(K));
    sM = reshape(ME(:) * denSld(:)', length(ME)*length(denSld), 1);
    M = fsparse(Iar(:,1), Iar(:,2), sM, [cfg.nDof, cfg.nDof]);
    M = M + M' - diag(diag(M));
    M(bc.massDof, bc.massDof) = M(bc.massDof, bc.massDof) + cfg.mass * eye(length(bc.massDof));
    M = M + fsparse(1:cfg.nDof, 1:cfg.nDof, eps*ones(1, cfg.nDof), [cfg.nDof, cfg.nDof]);
    K = K + fsparse(1:cfg.nDof, 1:cfg.nDof, eps*ones(1, cfg.nDof), [cfg.nDof, cfg.nDof]);

    [PsiAll, LambAll] = eigs_wrapper(K, M, freedofLft, cfg.n_ord_all);
    Lamb1 = LambAll(cfg.n_ord, cfg.n_ord);
    Psi1 = zeros(cfg.nDof,1); Psi1(freedofLft) = PsiAll(:, cfg.n_ord);
    Psi1 = Psi1 / sqrt(Psi1' * M * Psi1);

    OBJ(iter) = Lamb1;
    fval = sum(den) * cfg.EW*cfg.EH*cfg.EL / (cfg.DW*cfg.DH*cfg.DL) - cfg.volfrac;
    CONS(iter) = fval + cfg.volfrac;

    % sensitivity
    [a1, a2] = calc_f0val(K, M, freedofLft, Psi1, cfg.nDof, Lamb1);
    delta_H = 3*(1-cfg.alpha)/(4*epsilon) * (1 - Phimax.^2/(epsilon^2));
    delta_H(abs(Phimax) > epsilon) = 0;
    energy = sum((a1(grid.edofMat) .* (KE - Lamb1 .* ME)) .* Psi1(grid.edofMat) - 0.5 * (a2 * Psi1' * ME) .* Psi1(grid.edofMat), 2);
    sEner = energy * ones(1, 8) / 8;
    engyNod = sparse(double(grid.eleNodesID(:)), 1, sEner(:));
    df0dx = zeros(1, nDsvb); dfdx = zeros(1, nDsvb);
    df0dx(actDsvb) = (engyNod .* delta_H)' * PhimaxDrvAct;
    dfdx(actDsvb)  = (grid.volNod .* delta_H)' * PhimaxDrvAct * cfg.EW*cfg.EH*cfg.EL / (cfg.DW*cfg.DH*cfg.DL);

    dgt = cfg.dgt0 - floor(log10([max(abs(df0dx(:))) max(abs(dfdx(:)))]));
    scl = max(abs(df0dx));
    f0val = Lamb1 / scl;
    df0dx = round(df0dx * 10.^dgt(1)) / 10.^dgt(1) / scl;
    dfdx  = round(dfdx  * 10.^dgt(2)) / 10.^dgt(2);

    [xmma, ~, ~, ~, ~, ~, ~, low, upp] = mmasub(cfg.m, nDsvb, iter, xval(:), xmin, xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, cfg.a0, cfg.a, cfg.c, cfg.d);
    xold2 = xold1; xold1 = xval; xval = xmma;

    if iter >= 5 && fval/cfg.volfrac < 1e-4
        objVr5 = abs(max(abs(OBJ(iter-4:iter) - mean(OBJ(iter-4:iter)))) / mean(OBJ(iter-4:iter)));
    end

    loop = [loop iter]; var_obj = [var_obj Lamb1];
    fprintf('It %4d  Obj %.4f  Vol %.4f  Ch %.4f  Time %.2f\n', iter, f0val*scl, fval, objVr5, etime(clock, Time_iteration));
    iter = iter + 1;
end

results.OBJ = OBJ; results.loop = loop; results.CONS = CONS;

end

function cfg = defaults3d()
    cfg.DL = 8; cfg.DW = 4; cfg.DH = 4;
    cfg.nelx = 80; cfg.nely = 40; cfg.nelz = 40;
    cfg.vInt = [2.5 0.2 0.2 0 asin(1/sqrt(6)) asin(1/sqrt(5))];
    cfg.volfrac = 0.2; cfg.E = 2e11; cfg.nu = 0.3; cfg.rho = 7800; cfg.mass = 1e6;
    cfg.n_ord = 1; cfg.n_ord_all = 6; cfg.p = 6; cfg.lmd = 100; cfg.maxiter = 500;
    cfg.alpha = 1e-6; cfg.mid_num = 30; cfg.rank = 1; cfg.ini = 0.1; cfg.ep_old = 0.2; cfg.dgt0 = 5;
    cfg.m = 1; cfg.c = 1000; cfg.d = 0; cfg.a0 = 1; cfg.a = 0;

    cfg.EL = cfg.DL / cfg.nelx; cfg.EW = cfg.DW / cfg.nely; cfg.EH = cfg.DH / cfg.nelz;
    cfg.minSz = min([cfg.EL, cfg.EW, cfg.EH]);

    % initial components
    cfg.x0 = [kron(cfg.DL/4:cfg.DL/2:cfg.DL, ones(1,8)), kron(cfg.DL/4:cfg.DL/2:cfg.DL, ones(1,8))];
    cfg.y0 = [repmat(kron(cfg.DW/4:cfg.DW/2:cfg.DW, ones(1,4)),1,2), repmat(kron(cfg.DW/4:cfg.DW/2:cfg.DW, ones(1,4)),1,2)];
    cfg.z0 = [cfg.DH/4*ones(1,16), 3*cfg.DH/4*ones(1,16)];
    cfg.N  = length(cfg.x0);
    cfg.l1 = repmat(cfg.vInt(1), 1, cfg.N);
    cfg.l2 = repmat(cfg.vInt(2), 1, cfg.N);
    cfg.l3 = repmat(cfg.vInt(3), 1, cfg.N);
    cfg.alp = repmat(cfg.vInt(4), 1, cfg.N);
    cfg.bet = repmat([1 -1 1 -1 1 -1 1 -1]*cfg.vInt(5), 1, cfg.N/8);
    cfg.gam = repmat([1 1 -1 -1 1 1 -1 -1]*cfg.vInt(6), 1, cfg.N/8);
    cfg.dd  = [cfg.x0 - cfg.DL/2; cfg.y0 - cfg.DW/2; cfg.z0 - cfg.DH/2; cfg.l1; cfg.l2; cfg.l3; cfg.alp; cfg.bet; cfg.gam];
    cfg.nNod = (cfg.nelx+1)*(cfg.nely+1)*(cfg.nelz+1);
    cfg.PhiNd = [];
    cfg.bcCase = 'appendix-g';

    cfg.xmin = [-cfg.DL/2; -cfg.DW/2; -cfg.DH/2; cfg.minSz*[1;1;1]; -pi; -pi; -pi];
    cfg.xmax = [ cfg.DL/2;  cfg.DW/2;  cfg.DH/2; sqrt(cfg.DL^2+cfg.DW^2+cfg.DH^2)/2*[1;1;1]; pi; pi; pi];
    cfg.nDof = 3 * cfg.nNod;
end

function cfg = override(cfg, user)
    f = fieldnames(user);
    for i = 1:numel(f)
        cfg.(f{i}) = user.(f{i});
    end
end

function selftest_3d()
%SELFTEST_3D quick 3-iteration sanity check
    cfg = struct('nelx',10,'nely',4,'nelz',4,'maxiter',3,'bcCase','appendix-g');
    mmc_nfreq_3d_main(cfg);
    fprintf('Selftest complete. Inspect eigenvalues for positivity and K/M symmetry.\n');
end


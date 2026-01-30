% MMC 3D natural frequency topology optimization (cleaned from Appendix G)
% OCR fixes applied (key corrections):
%  - (nelx+1i) -> (nelx+1); nodMat reshape dims fixed to (nely+1,nelx+1,nelz+1)
%  - edofMat assembly restored: edofVec + [3*nNodfc+[4 5 6 1 2 3] -2 -1 0 1 2 3]
%  - massEle, fixNd/fixEle indexing corrected per appendix; stray minus signs removed
%  - variables ik/jK, Iar/Iar0, etc. standardized; missing ends/brackets added
%  - reshape sizes corrected (e.g., KE 24x24, ME 24x24), calc_Phi temp expression fixed

clear; clc; close all;
fsparse = @(i,j,s,sz) sparse(i,j,s,sz(1),sz(2));

% --------------------------- SEC 1): PARAMETERS SETTING
DL = 8; DW = 4; DH = 4;
nelx = 80; nely = 40; nelz = 40;
vInt = [2.5 0.2 0.2 0 asin(1/sqrt(6)) asin(1/sqrt(5))];
volfrac = 0.2;
E = 2e11; nu = 0.3; rho = 7800;
mass = 1e6;
n_ord = 1; n_ord_all = 6;
dgt0 = 5; scl_old = 1;
p = 6; lmd = 100;
iter = 1; maxiter = 500;
objVr5 = 1.0;
alpha = 1e-6; mid_num = 30; rank = 1; ini = 0.1; ep_old = 0.2;
loop = []; var_obj = []; Lamb_w1 = []; Lamb_w2 = [];

% --------------------------- SEC 2): SETTING OF FE DISCRETIZATION
nEle = nelx*nely*nelz;
nNod = (nelx+1)*(nely+1)*(nelz+1);
nNodfc = (nelx+1)*(nely+1);
nDof = 3*nNod;
EL = DL/nelx; EW = DW/nely; EH = DH/nelz;
minSz = min([EL, EW, EH]);
[Ke] = Ke_tril(E, nu, EL, EW, EH);
KE = zeros(24);
KE(tril(ones(24)) == 1) = Ke';
KE = KE + KE' - diag(diag(KE));
[Me] = Me_tril(rho, EL, EW, EH);
ME = zeros(24);
ME(tril(ones(24)) == 1) = Me';
ME = ME + ME' - diag(diag(ME));
ltMask = tril(true(24));
KElt = KE(ltMask);
MElt = ME(ltMask);
nodMat = int32(reshape(1:nNod, nely+1, nelx+1, nelz+1));
edofVec = reshape(3*nodMat(1:nely, 1:nelx, 1:nelz), nEle, 1);
base = [0, 3, 3*((nely+1)+1), 3*(nely+1), ...
        3*nNodfc, 3*(nNodfc+1), 3*(nNodfc+(nely+1)+1), 3*(nNodfc+(nely+1))];
offsets = int32(kron(base, ones(1,3)) + kron(ones(1,8), [-2 -1 0]));
edofMat = edofVec + offsets;
eleNodesID = edofMat(:, 3:3:24) ./ 3;
[sI, sII] = deal([]);
for j = 1:24
    sI  = cat(2, sI, j:24);
    sII = cat(2, sII, repmat(j, 1, 24-j+1));
end
[iK, jK] = deal(edofMat(:, sI)', edofMat(:, sII)');
Iar0 = sort([iK(:), jK(:)], 2, 'descend'); clear iK jK
[x, y, z] = meshgrid(EL*(-nelx/2:nelx/2), EW*(-nely/2:nely/2), EH*(-nelz/2:nelz/2));
LSgrid.x = permute(x, [2,1,3]); LSgrid.y = permute(y, [2,1,3]); LSgrid.z = permute(z, [2,1,3]);
volNod = sparse(double(eleNodesID(:)), 1, 1/8);

% --------------------------- SEC 3): LOADS, DISPLACEMENT BOUNDARY CONDITIONS
[jN, kN] = meshgrid(1:nely+1, 1:nelz+1);
fixNd = [1; (nelx+1)*(nely+1); nelx+1; nely*(nelx+1)+1];
fixDof = [3*fixNd(:)-2; 3*fixNd(:)-1; 3*fixNd(:)];
[jE, kE] = meshgrid(1:nely, 1:nelz);
fixEle = [1; nelx; nelx*nely; nelx*(nely-1)+1; (nelz-1)*nelx*nely+1; (nelz-1)*nelx*nely+nelx; (nelz-1)*nelx*nely+nelx*nely; (nelz-1)*nelx*nely+nelx*(nely-1)+1];
freeDof = setdiff(1:nDof, fixDof);
massNd = (nely/2)*(nelx+1) + (nelx/2) + 1 + nelz*((nelx+1)*(nely+1));
massDof = 3*massNd-2:3*massNd;
massEle = (nelx/2) + (nely/2 - 1)*nelx + (nelz - 1)*nelx*nely;

% --------------------------- SEC 4): INITIAL SETTING OF COMPONENTS
x0 = [kron(DL/4:DL/2:DL, ones(1,8)), kron(DL/4:DL/2:DL, ones(1,8))];
y0 = [repmat(kron(DW/4:DW/2:DW, ones(1,4)), 1, 2), repmat(kron(DW/4:DW/2:DW, ones(1,4)), 1, 2)];
z0 = [DH/4*ones(1,16), 3*DH/4*ones(1,16)];
N  = length(x0);
l1 = repmat(vInt(1), 1, N);
l2 = repmat(vInt(2), 1, N);
l3 = repmat(vInt(3), 1, N);
alp = repmat(vInt(4), 1, N);
bet = repmat([1 -1 1 -1 1 -1 1 -1]*vInt(5), 1, N/8);
gam = repmat([1 1 -1 -1 1 1 -1 -1]*vInt(6), 1, N/8);
dd  = [x0-DL/2; y0-DW/2; z0-DH/2; l1; l2; l3; alp; bet; gam];
nDsvb = length(dd(:));
nEhcp = nDsvb/N;
actComp = 1:N;
actDsvb = 1:nDsvb;
nNd = 0; PhiNd = [];
allPhi = [zeros(nNod, N) PhiNd];

% --------------------------- SEC 5): SETTING OF MMA
m = 1; c = 1000*ones(m,1); d = zeros(m,1);
a0 = 1; a = zeros(m,1);
xval = dd(:); xold1 = xval; xold2 = xval;
xmin = [-DL/2; -DW/2; -DH/2; minSz*[1;1;1]; -pi; -pi; -pi];
xmax = [ DL/2;  DW/2;  DH/2; sqrt(DL^2+DW^2+DH^2)/2*[1;1;1];  pi;  pi;  pi];
xmin = repmat(xmin, N, 1); xmax = repmat(xmax, N, 1);
low = xmin; upp = xmax;

% --------------------------- SEC 6): OPTIMIZATION LOOP
while objVr5 > 1e-4 && iter <= maxiter
    Time_iteration = clock;
    epsilon = ini + ep_old - (1/ep_old + exp(-rank*(iter - mid_num)))^(-1);

    % LP1: Generating TDFs and their derivatives
    allPhiDrv = sparse(nNod, nDsvb);
    for i = actComp
        [allPhi, allPhiDrv, xval, actComp, actDsvb] = ...
            calc_Phi(allPhi, allPhiDrv, xval, i, LSgrid, p, nEhcp, epsilon, actComp, actDsvb, minSz);
    end
    allPhiAct = [allPhi(:, actComp) PhiNd];
    temp = exp(lmd * allPhiAct);
    Phimax = max(-1e3, log(sum(temp, 2)) / lmd);
    allPhiDrvAct = allPhiDrv(:, actDsvb);
    Phimaxdphi = kron(temp(:, 1:length(actComp)) ./ (sum(temp, 2) + eps), ones(1, nEhcp));
    PhimaxDrvAct = Phimaxdphi .* allPhiDrvAct;

    % LP2: Plotting current design (optional)
    % figure(1); clf;
    % h = patch(isosurface(x, y, z, permute(reshape(Phimax, nelx+1, nely+1, nelz+1), [2,1,3]), 0));
    % hi = patch(isocaps(x, y, z, permute(reshape(Phimax, nelx+1, nely+1, nelz+1), [2,1,3]), 0));
    % set(h, 'FaceColor', 'red', 'EdgeColor', 'none'); set(hi, 'FaceColor','interp','EdgeColor','none');
    % colormap([1 0 0]);
    % isonormals(x, y, z, permute(reshape(Phimax, nelx+1, nely+1, nelz+1), [2,1,3]), h);
    % lighting flat; view(3); axis image; axis([-DL/2, DL/2, -DW/2, DW/2, -DH/2, DH/2]); light; drawnow;

    % LP3: Finite element analysis
    H = Heaviside(Phimax, alpha, epsilon);
    den = sum(H(eleNodesID), 2) / 8;
    Psi1 = zeros(nDof, 1); Psi2 = zeros(nDof, 1);
    bw_den = reshape(den, nelx, nely, nelz) > alpha; bw_den_sol = reshape(bwlabeln(bw_den), nelx*nely*nelz, 1);
    bw_fix = unique(bw_den_sol(fixEle)); bw_mass = nonzeros(unique(bw_den_sol(massEle)));

    struct = 0; denSld = den; edofMatLft = edofMat; freedofLft = freeDof; Iar = Iar0;
    if sum(ismember(bw_fix, bw_mass)) > 0
        struct = 1;
        eleLft = find(bw_den_sol == bw_mass);
        denSld = zeros(size(den)); denSld(eleLft) = den(eleLft);
        edofMatLft = edofMat(eleLft, :);
        freedofLft = setdiff(edofMatLft(:), fixDof);
        [iK1, jK1] = deal(edofMatLft(:, sI)', edofMatLft(:, sII)');
        Iar = sort([iK1(:), jK1(:)], 2, 'descend'); clear iK1 jK1
    end

    sK = kron(denSld(:), KElt(:));
    K = fsparse(Iar(:,1), Iar(:,2), sK, [nDof, nDof]); K = K + K' - diag(diag(K));
    sM = kron(denSld(:), MElt(:));
    M = fsparse(Iar(:,1), Iar(:,2), sM, [nDof, nDof]); M = M + M' - diag(diag(M));
    M(massDof, massDof) = M(massDof, massDof) + mass * eye(length(massDof));
    M = M + fsparse(1:nDof, 1:nDof, eps*ones(1, nDof), [nDof, nDof]);
    K = K + fsparse(1:nDof, 1:nDof, eps*ones(1, nDof), [nDof, nDof]);

    Time_FEA = clock;
    Kloc = (K(freedofLft, freedofLft) + K(freedofLft, freedofLft)')/2;
    Mloc = (M(freedofLft, freedofLft) + M(freedofLft, freedofLft)')/2;
    try
        [PsiAll, LambAll] = eigs(Kloc, Mloc, n_ord_all, 'sm');
    catch
        reg = 1e-6;
        [PsiAll, LambAll] = eigs(Kloc + reg*speye(size(Kloc,1)), Mloc, n_ord_all, 'sm');
    end
    Lamb1 = LambAll(n_ord, n_ord);
    Psi1(freedofLft) = PsiAll(:, n_ord);
    Psi1 = Psi1 / sqrt(Psi1' * M * Psi1);

    OBJ(iter) = Lamb1;
    fval = sum(den) * EW * EH * EL / (DW * DH * DL) - volfrac; CONS(iter) = fval + volfrac;

    % LP4: Sensitivity analysis
    df0dx = zeros(1, nDsvb); dfdx = zeros(1, nDsvb);
    delta_H = 3*(1-alpha)/(4*epsilon) * (1 - Phimax.^2 / (epsilon^2));
    delta_H(abs(Phimax) > epsilon) = 0;
    [a1, a2] = calc_f0val(K, M, freedofLft, Psi1, nDof, Lamb1);
    a1e  = a1(edofMat);
    psiE = Psi1(edofMat);
    term1 = (a1e * (KE - Lamb1*ME)) .* psiE;
    term2 = 0.5 * a2 * (psiE * ME) .* psiE;
    energy = sum(term1 - term2, 2);
    sEner = energy * ones(1, 8) / 8;
    engyNod = sparse(double(eleNodesID(:)), 1, sEner(:));
    df0dx(actDsvb) = (engyNod .* delta_H)' * PhimaxDrvAct;
    dfdx(actDsvb) = (volNod .* delta_H)' * PhimaxDrvAct * EW * EH * EL / (DW * DH * DL);

    dgt = dgt0 - floor(log10([max(abs(df0dx(:))) max(abs(dfdx(:)))]));
    scl = max(abs(df0dx));
    f0val = Lamb1 / scl;
    df0dx = round(df0dx * 10.^dgt(1)) / 10.^dgt(1) / scl;
    dfdx  = round(dfdx  * 10.^dgt(2)) / 10.^dgt(2);

    df0dx_col = df0dx(:);
    dfdx_mat  = reshape(dfdx, 1, nDsvb);
    fval_vec  = fval(:);
    [xmma, ~, ~, ~, ~, ~, ~, low, upp] = mmasub(m, nDsvb, iter, xval(:), xmin, xmax, xold1, xold2, f0val, df0dx_col, fval_vec, dfdx_mat, low, upp, a0, a, c, d);
    xold2 = xold1; xold1 = xval; xval = xmma;

    if iter >= 5 && fval/volfrac < 1e-4
        objVr5 = abs(max(abs(OBJ(iter-4:iter) - mean(OBJ(iter-4:iter)))) / mean(OBJ(iter-4:iter)));
    end

    loop = [loop iter]; var_obj = [var_obj Lamb1];
    fprintf('? It.: %4i\t Obj.: %6.3f\t Vol.: %6.4f\t Ch.: %6.4f\t Time_iteration: %6.4f\n', iter, f0val*scl, fval, objVr5, etime(clock, Time_iteration));
    iter = iter + 1;
end

% --------------------------- functions ---------------------------
function [a1, a2] = calc_f0val(K, M, Dof, Psi, nDof, Lamb)
L11 = K(Dof, Dof) - Lamb * M(Dof, Dof);
L12 = -M(Dof, Dof) * Psi(Dof);
rhs = [zeros(length(Dof),1); 1];
LL  = [L11, L12; L12', 0];
vec_a = LL \ rhs;
a1 = zeros(nDof, 1); a1(Dof) = vec_a(1:length(Dof));
a2 = vec_a(end);
end

function [allPhi, allPhidrv, xval, actComp, actDsvb] = calc_Phi(allPhi, allPhidrv, xval, i, LSgrid, p, nEhcp, epsilon, actComp, actDsvb, minSz)
di = xval((i-1)*nEhcp+1 : i*nEhcp);
x0 = di(1); y0 = di(2); z0 = di(3); l1 = di(4) + eps; l2 = di(5) + eps; l3 = di(6) + eps;
sa = sin(di(7)); sb = sin(di(8)); sg = sin(di(9));
ca = cos(di(7)); cb = cos(di(8)); cg = cos(di(9));
R = [cb*cg, cb*sg, -sb;
     sa*sb*cg - ca*sg, sa*sb*sg + ca*cg, sa*cb;
     ca*sb*cg + sa*sg, ca*sb*sg - sa*cg, ca*cb];
xyzLc = [LSgrid.x(:)-x0, LSgrid.y(:)-y0, LSgrid.z(:)-z0];
xyz   = xyzLc * R';
x1 = xyz(:,1) + eps; y1 = xyz(:,2) + eps; z1 = xyz(:,3) + eps;
temp = (abs(x1./l1).^p + abs(y1./l2).^p + abs(z1./l3).^p);
allPhi(:, i) = 1 - temp.^(1/p);

if (l1/minSz < 1.01 && l2/minSz < 1.01) || min(abs(allPhi(:, i))) >= epsilon
    allPhi(:, i) = -1e3;
    xval((i-1)*nEhcp + (4:6)) = 0;
    actComp = setdiff(actComp, i);
    actDsvb = setdiff(actDsvb, nEhcp*i - nEhcp + 1 : nEhcp*i);
    return;
end

Ra = [0 0 0; R(3,1) R(3,2) R(3,3); -R(2,1) -R(2,2) -R(2,3)];
Rb = [-sb*cg -sb*sg -cb; sa*cb*cg sa*cb*sg -sa*sb; ca*cb*cg ca*cb*sg -ca*sb];
Rg = [-cb*sg cb*cg 0; -sa*sb*sg - ca*cg sa*sb*cg - sa*sg 0; -ca*sb*sg + sa*cg ca*sb*cg + sa*sg 0];

dxi = xyzLc * [Ra(1,:); Rb(1,:); Rg(1,:)];
dyi = xyzLc * [Ra(2,:); Rb(2,:); Rg(2,:)];
dzi = xyzLc * [Ra(3,:); Rb(3,:); Rg(3,:)];

temp1 = -temp.^(1/p - 1) .* (x1./l1).^(p-1) / l1;
temp2 = -temp.^(1/p - 1) .* (y1./l2).^(p-1) / l2;
temp3 = -temp.^(1/p - 1) .* (z1./l3).^(p-1) / l3;

dpdxi = temp1; dpdyi = temp2; dpdzi = temp3;
dpdl1 = -temp1 .* (x1./l1); dpdl2 = -temp2 .* (y1./l2); dpdl3 = -temp3 .* (z1./l3);

% Derivatives w.r.t design vars [x0 y0 z0 l1 l2 l3 alpha beta gamma]
Phi_drv_block = [ -dpdxi, -dpdyi, -dpdzi, dpdl1, dpdl2, dpdl3, zeros(length(x1),3) ];
allPhidrv(:, nEhcp*(i-1)+1 : nEhcp*i) = Phi_drv_block;
end

function H = Heaviside(phi, alpha, epsilon)
H = 3*(1-alpha)/4 * (phi/epsilon - phi.^3/(3*(epsilon)^3)) + (1+alpha)/2;
H(phi > epsilon)  = 1;
H(phi < -epsilon) = alpha;
end

function Ke = Ke_tril(E, nu, a, b, h)
% Wrapper: use validated HEX8 formulation from fe_hex8, return lower tri vector
    [KEfull, ~] = fe_hex8(E, nu, 1, a, b, h);  % rho ignored here
    Ke = KEfull(tril(ones(24)) == 1)';
end

function Me = Me_tril(rho, a, b, h)
% Wrapper: use validated HEX8 mass, return lower tri vector
    [~, MEfull] = fe_hex8(1, 0, rho, a, b, h);  % E,nu dummy
    Me = MEfull(tril(ones(24)) == 1)';
end

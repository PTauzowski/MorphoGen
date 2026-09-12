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
nelx = 40; nely = 20; nelz = 20;
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

% Debug: check element matrices
fprintf('DEBUG: Element dimensions: EL=%.3f, EW=%.3f, EH=%.3f m\n', EL, EW, EH);
fprintf('DEBUG: Element volume = %.6f m^3\n', EL*EW*EH);
fprintf('DEBUG: max(KE) = %.3e, trace(KE) = %.3e\n', max(KE(:)), trace(KE));
fprintf('DEBUG: max(ME) = %.3e, trace(ME) = %.3e\n', max(ME(:)), trace(ME));
fprintf('DEBUG: Element mass from ME (trace/3) = %.3f kg\n', trace(ME)/3);  % Approx
fprintf('DEBUG: Expected element mass (rho*V) = %.3f kg\n', rho*EL*EW*EH);
ltMask = tril(true(24));
KElt = KE(ltMask);
MElt = ME(ltMask);
Ke = KElt;  % Lower triangular vector for assembly
Me = MElt;  % Lower triangular vector for assembly
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
% Fixed nodes at 4 corners of z=0 face (bottom layer)
% Node numbering: nodMat(i,j,k) where i=y-index, j=x-index, k=z-index
% Linear index = i + (j-1)*(nely+1) + (k-1)*(nely+1)*(nelx+1)
% Corner (0,0,0): nodMat(1,1,1) = 1
% Corner (DL,0,0): nodMat(1,nelx+1,1) = 1 + nelx*(nely+1)
% Corner (0,DW,0): nodMat(nely+1,1,1) = nely+1
% Corner (DL,DW,0): nodMat(nely+1,nelx+1,1) = (nelx+1)*(nely+1)
fixNd = [1; nely+1; 1+nelx*(nely+1); (nelx+1)*(nely+1)];
fixDof = [3*fixNd(:)-2; 3*fixNd(:)-1; 3*fixNd(:)];  % Fix all 3 DOFs at corners
freeDof = setdiff(1:nDof, fixDof);

% Fixed elements at 4 corners of z=0 face
% Element numbering: ele = i + (j-1)*nely + (k-1)*nely*nelx
% Corner elements at z=0: (1,1), (nely,1), (1,nelx), (nely,nelx)
fixEle_z0 = [1; nely; (nelx-1)*nely+1; nelx*nely];
fixEle = fixEle_z0;  % Only z=0 corners (supports are at bottom)

% Mass at center of top face (z=DH)
% massNd at (x=DL/2, y=DW/2, z=DH) -> node indices (nely/2+1, nelx/2+1, nelz+1)
massNd = (nely/2+1) + (nelx/2)*(nely+1) + nelz*(nely+1)*(nelx+1);
massDof = 3*massNd-2:3*massNd;

% massEle: element containing mass node (at top layer, center)
% Element at (i=nely/2, j=nelx/2, k=nelz)
massEle = (nely/2) + (nelx/2-1)*nely + (nelz-1)*nely*nelx;

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

    % LP2: Plotting current design
    figure(1); clf;
    PhiPlot = permute(reshape(Phimax, nelx+1, nely+1, nelz+1), [2,1,3]);
    h = patch(isosurface(x, y, z, PhiPlot, 0));
    hi = patch(isocaps(x, y, z, PhiPlot, 0));
    set(h, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    set(hi, 'FaceColor', 'interp', 'EdgeColor', 'none');
    colormap([0.8 0.2 0.2]);
    isonormals(x, y, z, PhiPlot, h);
    lighting gouraud; view(3); axis equal;
    axis([-DL/2, DL/2, -DW/2, DW/2, -DH/2, DH/2]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    light('Position', [1 1 1]); light('Position', [-1 -1 0.5]);
    title(sprintf('Iteration %d', iter));
    drawnow;

    % LP3: Finite element analysis
    H = Heaviside(Phimax, alpha, epsilon);
    den = sum(H(eleNodesID), 2) / 8;
    Psi1 = zeros(nDof, 1); Psi2 = zeros(nDof, 1);
    % Note: element ordering follows nodMat indexing (nely, nelx, nelz)
    bw_den = reshape(den, nely, nelx, nelz) > alpha;
    bw_den_sol = reshape(bwlabeln(bw_den), nely*nelx*nelz, 1);
    bw_fix = unique(bw_den_sol(fixEle));
    bw_mass = nonzeros(unique(bw_den_sol(massEle)));

    % Debug output for first iteration
    if iter == 1
        fprintf('DEBUG: fixEle = [%s]\n', num2str(fixEle'));
        fprintf('DEBUG: massEle = %d\n', massEle);
        fprintf('DEBUG: bw_fix labels = [%s]\n', num2str(bw_fix'));
        fprintf('DEBUG: bw_mass labels = [%s]\n', num2str(bw_mass'));
        fprintf('DEBUG: Number of connected components = %d\n', max(bw_den_sol));
        fprintf('DEBUG: Volume of solid elements = %.4f\n', sum(bw_den_sol > 0)/nEle);
    end

    clear bw_den
    struct = 0;
    if sum(ismember(bw_fix, bw_mass)) > 0
        struct = 1;
        eleLft = find(bw_den_sol == bw_mass);
        denSld = zeros(size(den)); denSld(eleLft) = den(eleLft);
        edofMatLft = edofMat(eleLft, :);
        freedofLft = setdiff(edofMatLft(:), fixDof);
        [iK1, jK1] = deal(edofMatLft(:, sI)', edofMatLft(:, sII)');
        Iar = sort([iK1(:), jK1(:)], 2, 'descend'); clear iK1 jK1

        sK = reshape(Ke(:)*denSld(eleLft)', length(Ke)*length(eleLft), 1);
        K = fsparse(Iar(:,1), Iar(:,2), sK, [nDof, nDof]); K = K + K' - diag(diag(K));
        K = K + fsparse(1:nDof, 1:nDof, eps*ones(1,nDof), [nDof, nDof]);
        sM = reshape(Me(:)*denSld(eleLft)', length(Me)*length(eleLft), 1);
        M = fsparse(Iar(:,1), Iar(:,2), sM, [nDof, nDof]); M = M + M' - diag(diag(M));
        M(massDof, massDof) = M(massDof, massDof) + mass*eye(length(massDof));
        M = M + fsparse(1:nDof, 1:nDof, eps*ones(1,nDof), [nDof, nDof]);
        Time_FEA = clock;
        Kred = K(freedofLft, freedofLft);
        Mred = M(freedofLft, freedofLft);

        % Debug: check matrix properties at iteration 1
        if iter == 1
            fprintf('DEBUG: Size of reduced K: %d x %d\n', size(Kred));
            fprintf('DEBUG: nnz(K) = %d, nnz(M) = %d\n', nnz(Kred), nnz(Mred));
            fprintf('DEBUG: max(K) = %.3e, min(K) = %.3e\n', full(max(Kred(:))), full(min(Kred(:))));
            fprintf('DEBUG: max(M) = %.3e, min(M) = %.3e\n', full(max(Mred(:))), full(min(Mred(:))));
        end

        [PsiAll, LambAll] = eigs(Kred, Mred, n_ord_all, 'sm');

        % Debug: show all computed eigenvalues at iteration 1
        if iter == 1
            eigvals = diag(LambAll);
            fprintf('DEBUG: Computed eigenvalues: ');
            fprintf('%.3e ', eigvals);
            fprintf('\n');
            fprintf('DEBUG: Corresponding frequencies (Hz): ');
            fprintf('%.2f ', sqrt(abs(eigvals))/(2*pi));
            fprintf('\n');
        end

        Lamb1 = LambAll(n_ord, n_ord);
        Psi1(freedofLft) = PsiAll(:, n_ord);
        Psi1 = Psi1 / sqrt(Psi1' * M * Psi1);
    else
        % NO loading path found - use regular FEA with all elements
        disp('WARNING!!! NO loading path is founded!!!');
        sK = reshape(Ke(:)*den(:)', length(Ke)*nEle, 1);
        K = fsparse(Iar0(:,1), Iar0(:,2), sK, [nDof, nDof]); K = K + K' - diag(diag(K));
        K = K + fsparse(1:nDof, 1:nDof, eps*ones(1,nDof), [nDof, nDof]);  % regularization
        sM = reshape(Me(:)*den(:)', length(Me)*nEle, 1);
        M = fsparse(Iar0(:,1), Iar0(:,2), sM, [nDof, nDof]); M = M + M' - diag(diag(M));
        M(massDof, massDof) = M(massDof, massDof) + mass*eye(length(massDof));
        M = M + fsparse(1:nDof, 1:nDof, eps*ones(1,nDof), [nDof, nDof]);
        Time_FEA = clock;
        try
            [PsiAll, LambAll] = eigs(K(freeDof, freeDof), M(freeDof, freeDof), n_ord_all, 'sm');
        catch
            % If eigs fails, add stronger regularization
            Kreg = K(freeDof, freeDof) + 1e-6*speye(length(freeDof));
            [PsiAll, LambAll] = eigs(Kreg, M(freeDof, freeDof), n_ord_all, 'sm');
        end
        Lamb1 = LambAll(n_ord, n_ord);
        Psi1(freeDof) = PsiAll(:, n_ord);
        Psi1 = Psi1 / sqrt(Psi1' * M * Psi1);
        freedofLft = freeDof;  % Use all free DOFs for sensitivity
    end

    OBJ(iter) = Lamb1;
    volCur = sum(den) * EW * EH * EL / (DW * DH * DL);  % Current volume fraction
    fval = volCur - volfrac;  % Constraint: vol - target (should be <= 0)
    CONS(iter) = volCur;

    % LP4: Sensitivity analysis
    df0dx = zeros(1, nDsvb); dfdx = zeros(1, nDsvb);
    delta_H = 3*(1-alpha)/(4*epsilon) * (1 - Phimax.^2 / (epsilon^2));
    delta_H(abs(Phimax) > epsilon) = 0;
    if struct == 1
        [a1_w1, a2_w1] = calc_f0val(K, M, freedofLft, Psi1, nDof, Lamb1);
    else
        [a1_w1, a2_w1] = calc_f0val(K, M, freeDof, Psi1, nDof, Lamb1);
    end
    energy = sum((a1_w1(edofMat)*(KE - Lamb1.*ME)).*Psi1(edofMat) - 0.5*(a2_w1*Psi1(edofMat)*ME).*Psi1(edofMat), 2);
    sEner = energy * ones(1, 8) / 8;
    engyNod = sparse(double(eleNodesID(:)), 1, sEner(:));
    df0dx(actDsvb) = (engyNod .* delta_H)' * PhimaxDrvAct;
    dfdx(actDsvb) = (volNod .* delta_H)' * PhimaxDrvAct * EW * EH * EL / (DW * DH * DL);

    dgt = dgt0 - floor(log10([max(abs(df0dx(:))) max(abs(dfdx(:)))]));
    % Adaptive scl criterion from paper (Eq. 24)
    if iter == 1
        scl = max(abs(df0dx)); scl_old = scl;
    elseif OBJ(iter) - OBJ(iter-1) > 0
        scl = max(abs(df0dx)); scl_old = scl;
    else
        scl = scl_old;
    end
    f0val = -Lamb1 / scl;  % Negative for maximization (MMA minimizes)
    df0dx = -round(df0dx * 10.^dgt(1)) / 10.^dgt(1) / scl;  % Negative for maximization
    dfdx  = round(dfdx  * 10.^dgt(2)) / 10.^dgt(2);

    df0dx_col = df0dx(:);
    dfdx_mat  = reshape(dfdx, 1, nDsvb);
    fval_vec  = fval(:);
    [xmma, ~, ~, ~, ~, ~, ~, low, upp] = mmasub(m, nDsvb, iter, xval(:), xmin, xmax, xold1, xold2, f0val, df0dx_col, fval_vec, dfdx_mat, low, upp, a0, a, c, d);
    xold2 = xold1; xold1 = xval; xval = xmma;

    if iter >= 5 && abs(fval)/volfrac < 1e-2
        objVr5 = abs(max(abs(OBJ(iter-4:iter) - mean(OBJ(iter-4:iter)))) / mean(OBJ(iter-4:iter)));
    end

    loop = [loop iter]; var_obj = [var_obj Lamb1];
    freqHz = sqrt(max(0, Lamb1)) / (2*pi);  % Natural frequency in Hz
    fprintf('It.: %4i\t Freq.: %8.2f Hz\t Vol.: %6.4f/%4.2f\t Conv.: %6.4f\t Time: %5.2f s\n', ...
        iter, freqHz, volCur, volfrac, objVr5, etime(clock, Time_iteration));
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

function [allPhi, allPhidrv, xval, actComp, actDsvb] = calc_Phi( ...
    allPhi, allPhidrv, xval, i, LSgrid, p, nEhcp, epsilon, actComp, actDsvb, minSz)
%CALC_PHI MMC level-set field for component i + derivatives.
%
% di = [x0 y0 z0 l1 l2 l3 alpha beta gamma]  (nEhcp = 9)
%
% IMPORTANT:
%   Do NOT deactivate components just because l ~= minSz.
%   Your MMA xmin already enforces l >= minSz.
%   Deactivate ONLY if component has zero influence everywhere:
%       max(phi) <= -epsilon
%
% Convention:
%   phi > 0 inside component, phi = 0 boundary, phi < 0 outside.

    % -------------------- unpack design variables --------------------
    di = xval((i-1)*nEhcp+1 : i*nEhcp);

    x0 = di(1);  y0 = di(2);  z0 = di(3);
    l1 = di(4) + eps;  l2 = di(5) + eps;  l3 = di(6) + eps;
    a  = di(7);  b  = di(8);  g  = di(9);

    sa = sin(a); ca = cos(a);
    sb = sin(b); cb = cos(b);
    sg = sin(g); cg = cos(g);

    % -------------------- rotation matrix (your convention) --------------------
    R = [ cb*cg,                 cb*sg,               -sb;
          sa*sb*cg - ca*sg,      sa*sb*sg + ca*cg,    sa*cb;
          ca*sb*cg + sa*sg,      ca*sb*sg - sa*cg,    ca*cb ];

    % -------------------- centered global coords --------------------
    xyzLc = [LSgrid.x(:) - x0, LSgrid.y(:) - y0, LSgrid.z(:) - z0];  % (N x 3)

    % rotated coords: xyz = xyzLc * R'
    xyz = xyzLc * R.';       % (N x 3)
    x1  = xyz(:,1) + eps;
    y1  = xyz(:,2) + eps;
    z1  = xyz(:,3) + eps;

    % -------------------- MMC super-ellipsoid field --------------------
    tp  = (abs(x1./l1).^p + abs(y1./l2).^p + abs(z1./l3).^p);
    phi = 1 - tp.^(1/p);

    allPhi(:, i) = phi;

    % -------------------- deactivation (robust) --------------------
    % Component has no influence if it is everywhere below -epsilon
    noInfluence = (max(phi) <= -epsilon);

    if noInfluence
        allPhi(:, i) = -1e3;  % remove from KS-max safely
        actComp = setdiff(actComp, i);
        dv_idx = nEhcp*(i-1)+1 : nEhcp*i;
        actDsvb = setdiff(actDsvb, dv_idx);
        return;
    end

    % -------------------- derivatives of phi --------------------
    % phi = 1 - tp^(1/p)
    % dphi/dx1 = -tp^(1/p-1) * |x1/l1|^(p-1) * sign(x1) / l1
    tp_pow = tp.^(1/p - 1);
    sx = sign(x1);  sy = sign(y1);  sz = sign(z1);

    dphi_dx1 = - tp_pow .* (abs(x1./l1).^(p-1)) .* sx ./ l1;
    dphi_dy1 = - tp_pow .* (abs(y1./l2).^(p-1)) .* sy ./ l2;
    dphi_dz1 = - tp_pow .* (abs(z1./l3).^(p-1)) .* sz ./ l3;

    % length derivatives
    dphi_dl1 = -dphi_dx1 .* (x1./l1);
    dphi_dl2 = -dphi_dy1 .* (y1./l2);
    dphi_dl3 = -dphi_dz1 .* (z1./l3);

    % -------------------- translation derivatives --------------------
    % xyz = (X - [x0 y0 z0]) * R' => dxyz/dx0 = -R(1,:), etc.
    dxyz_dx0 = -R(1,:);
    dxyz_dy0 = -R(2,:);
    dxyz_dz0 = -R(3,:);

    dphi_dx0 = dphi_dx1*dxyz_dx0(1) + dphi_dy1*dxyz_dx0(2) + dphi_dz1*dxyz_dx0(3);
    dphi_dy0 = dphi_dx1*dxyz_dy0(1) + dphi_dy1*dxyz_dy0(2) + dphi_dz1*dxyz_dy0(3);
    dphi_dz0 = dphi_dx1*dxyz_dz0(1) + dphi_dy1*dxyz_dz0(2) + dphi_dz1*dxyz_dz0(3);

    % -------------------- angle derivatives --------------------
    % xyz = xyzLc * R' => dxyz/dtheta = xyzLc * (dR/dtheta)'
    dRda = [ 0, 0, 0;
             ca*sb*cg + sa*sg,   ca*sb*sg - sa*cg,    ca*cb;
            -sa*sb*cg + ca*sg,  -sa*sb*sg - ca*cg,   -sa*cb ];

    dRdb = [ -sb*cg,            -sb*sg,              -cb;
              sa*cb*cg,          sa*cb*sg,           -sa*sb;
              ca*cb*cg,          ca*cb*sg,           -ca*sb ];

    dRdg = [ -cb*sg,             cb*cg,               0;
             -sa*sb*sg - ca*cg,  sa*sb*cg - ca*sg,    0;
             -ca*sb*sg + sa*cg,  ca*sb*cg + sa*sg,    0 ];

    dxyz_da = xyzLc * dRda.';
    dxyz_db = xyzLc * dRdb.';
    dxyz_dg = xyzLc * dRdg.';

    dphi_da = dphi_dx1.*dxyz_da(:,1) + dphi_dy1.*dxyz_da(:,2) + dphi_dz1.*dxyz_da(:,3);
    dphi_db = dphi_dx1.*dxyz_db(:,1) + dphi_dy1.*dxyz_db(:,2) + dphi_dz1.*dxyz_db(:,3);
    dphi_dg = dphi_dx1.*dxyz_dg(:,1) + dphi_dy1.*dxyz_dg(:,2) + dphi_dz1.*dxyz_dg(:,3);

    % -------------------- write derivative block --------------------
    % Order: [x0 y0 z0 l1 l2 l3 alpha beta gamma]
    allPhidrv(:, nEhcp*(i-1)+1 : nEhcp*i) = ...
        [dphi_dx0, dphi_dy0, dphi_dz0, dphi_dl1, dphi_dl2, dphi_dl3, dphi_da, dphi_db, dphi_dg];

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

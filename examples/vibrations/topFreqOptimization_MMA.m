function [omega_best, xPhys_best] = topFreqOptimization_MMA( ...
    L, H, nelx, nely, volfrac, penal, rmin, maxiter, supportType)

% MMA topology optimization: maximize fundamental eigenfrequency
% (equivalently maximize lambda1 of K phi = lambda M phi).
% Reports omega = sqrt(lambda1) [rad/s].
%
% Uses:
% - density filter (radius rmin in element units)
% - Heaviside projection with beta continuation
% - MMA with one constraint: mean(xPhys) <= volfrac
%
% Requires mmasub.m and subsolv.m on MATLAB path.

%% -------------------- Material (paper-like) --------------------
E0      = 1.0e7;
Emin    = 1.0e-2;
rho0    = 1.0;
rho_min = 1.0e-5;
nu      = 0.3;
t       = 1.0;   % thickness

%% -------------------- Geometry / mesh -------------------------
dx = L/nelx;
dy = H/nely;
nEl = nelx*nely;

nodeNrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
nDof    = 2*(nelx+1)*(nely+1);

fixedDOFs = buildSupports(string(supportType), nodeNrs);
freeDOFs  = setdiff(1:nDof, fixedDOFs);

%% -------------------- Element matrices (Q4 rectangle) ----------
[Ke0, M0] = q4_rect_KeM_planeStress(E0, nu, rho0, t, dx, dy);
Ke_unit = Ke0 / E0;
M_unit  = M0  / rho0;

%% -------------------- Filter kernel ----------------------------
[dyf, dxf] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, -ceil(rmin)+1:ceil(rmin)-1);
h  = max(0, rmin - sqrt(dxf.^2 + dyf.^2));
bc = 'symmetric';
Hs = imfilter(ones(nely, nelx), h, bc);

% Helper: forward filter and "transpose" filter (same for symmetric kernel)
filtFwd  = @(xmat) imfilter(xmat, h, bc) ./ Hs;
filtBack = @(gmat) imfilter(gmat ./ Hs, h, bc);

%% -------------------- Projection (Heaviside) -------------------
eta     = 0.5;
betaMax = 16;

%% -------------------- Assembly indices -------------------------
cVec = 2*nodeNrs(1:nely,1:nelx) + 1;
cVec = cVec(:);
cMat = [cVec, cVec+1, ...
        cVec + 2*nely + 2, cVec + 2*nely + 3, ...
        cVec + 2*nely,     cVec + 2*nely + 1, ...
        cVec - 2,          cVec - 1];
cMat = double(cMat);

[Iloc,Jloc] = find(tril(ones(8)));
K0_lower = Ke_unit(sub2ind([8,8], Iloc, Jloc));
M0_lower = M_unit(sub2ind([8,8], Iloc, Jloc));

iK = reshape(cMat(:,Iloc)', [], 1);
jK = reshape(cMat(:,Jloc)', [], 1);

%% -------------------- MMA setup --------------------------------
n = nEl;    % number of design variables
m = 1;      % number of constraints (volume)

xmin = zeros(n,1);
xmax = ones(n,1);

xval  = volfrac*ones(n,1);
xold1 = xval;
xold2 = xval;

low = xmin;
upp = xmax;

a0 = 1;
a  = zeros(m,1);
c  = 1000*ones(m,1);
d  = zeros(m,1);

%% -------------------- Loop -------------------------------------
omega_best = -Inf;
xPhys_best = xval;

for iter = 1:maxiter

    % beta continuation: 1,2,4,8,16 every 25 iters
    beta = min(betaMax, 2^(floor((iter-1)/25)));

    % --- Filter + projection
    xTildeMat = filtFwd(reshape(xval,nely,nelx));
    [xPhysMat, dH] = heavisideProjection(xTildeMat, beta, eta);
    xPhys = xPhysMat(:);

    % --- Assemble K,M
    Ee = Emin + (xPhys.^penal) * (E0 - Emin);
    re = rho_min + xPhys * (rho0 - rho_min);

    Kval = kron(Ee, K0_lower);
    Mval = kron(re, M0_lower);

    K = sparse(iK, jK, Kval, nDof, nDof); K = K + K' - diag(diag(K));
    M = sparse(iK, jK, Mval, nDof, nDof); M = M + M' - diag(diag(M));

    Kf = K(freeDOFs,freeDOFs);
    Mf = M(freeDOFs,freeDOFs);

    % --- Eigen solve (smallest eigenvalue)
    opts.tol = 1e-9;
    opts.maxit = 800;
    opts.disp = 0;

    [vf, df] = eigs(Kf, Mf, 1, 'SM', opts);
    lambda1 = real(df(1,1));
    lambda1 = max(lambda1, 0);
    omega1  = sqrt(lambda1);

    if omega1 > omega_best
        omega_best = omega1;
        xPhys_best = xPhys;
    end

    % --- Mass normalize mode
    vf = vf / sqrt(real(vf'*(Mf*vf)));
    phi = zeros(nDof,1); phi(freeDOFs) = vf;

    % --- Sensitivities wrt xPhys: d lambda / d xPhys
    phi_e = phi(cMat);                        % nEl x 8
    phiK  = sum((phi_e * Ke_unit).*phi_e, 2); % phi^T Kunit phi
    phiM  = sum((phi_e * M_unit ).*phi_e, 2); % phi^T Munit phi

    dlam_dxPhys = (penal*(E0-Emin) * (xPhys.^(penal-1))) .* phiK ...
                - lambda1 * (rho0-rho_min) .* phiM;

    % Objective for MMA is MINIMIZATION:
    % f0 = -lambda  => maximize lambda (and omega)
    f0val = -lambda1;
    df0_dxPhys = -dlam_dxPhys;

    % Chain back: x -> filter -> projection -> xPhys
    df0_mat = reshape(df0_dxPhys, nely, nelx) .* dH;
    df0 = filtBack(df0_mat);
    df0 = df0(:);

    % Constraint: g(x) = mean(xPhys) - volfrac <= 0
    gval = mean(xPhys) - volfrac;
    dg_dxPhys = (1/nEl) * ones(nEl,1);

    dg_mat = reshape(dg_dxPhys, nely, nelx) .* dH;
    dg = filtBack(dg_mat);
    dg = dg(:);

    fval  = gval;             % m x 1
    dfdx  = dg';              % m x n  (MMA expects row(s))

    % --- MMA update
    [xnew,~,~,~,~,~,~,~,low,upp] = mmasub( ...
        m,n,iter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0,fval,dfdx,low,upp,a0,a,c,d);

    xold2 = xold1;
    xold1 = xval;
    xval  = xnew;

    % --- Print + plot
    fprintf('It:%3d  beta:%2d  omega: %.4f  best: %.4f  Vol(xPhys): %.3f  g: %+ .3e  BC:%s\n', ...
        iter, beta, omega1, omega_best, mean(xPhys), gval, upper(string(supportType)));

    colormap(gray);
    imagesc(1-reshape(xPhys,nely,nelx)); caxis([0 1]); axis equal off; drawnow;
end

fprintf('Best omega = %.4f rad/s\n', omega_best);

end

% -------------------------------------------------------------------------
function [xPhys, dH] = heavisideProjection(xTilde, beta, eta)
    denom = tanh(beta*eta) + tanh(beta*(1-eta));
    xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / denom;
    dH    = (beta * (1 - tanh(beta*(xTilde-eta)).^2)) / denom;
end

% -------------------------------------------------------------------------
function fixedDOFs = buildSupports(supportType, nodeNrs)
    nely = size(nodeNrs,1)-1;

    leftNodes  = nodeNrs(:,1);
    rightNodes = nodeNrs(:,end);

    leftBot  = nodeNrs(1,1);
    rightBot = nodeNrs(1,end);

    midRow   = floor((nely+1)/2) + 1;
    leftMid  = nodeNrs(midRow,1);
    rightMid = nodeNrs(midRow,end);

    u = @(n) 2*n - 1;
    v = @(n) 2*n;

    s = upper(string(supportType));

    switch s
        case "CF"
            fixedDOFs = [u(leftNodes(:)); v(leftNodes(:))];

        case "CC"
            fixedDOFs = [u(leftNodes(:)); v(leftNodes(:)); u(rightNodes(:)); v(rightNodes(:))];

        case "CH"
            fixedDOFs = [u(leftNodes(:)); v(leftNodes(:)); u(rightBot); v(rightBot)];

        case "HH"
            % Typical "pin + roller" model (simply supported):
            % left bottom: u=v=0, right bottom: v=0
            fixedDOFs = [u(leftBot); v(leftBot); v(rightBot)];

        case "HHMID"
            fixedDOFs = [u(leftMid); v(leftMid); v(rightMid)];

        otherwise
            error("Unknown supportType '%s'. Use CF, CC, CH, HH, HHmid.", s);
    end

    fixedDOFs = unique(fixedDOFs(:));
end

% -------------------------------------------------------------------------
function [Ke, Me] = q4_rect_KeM_planeStress(E, nu, rho, t, dx, dy)
    C = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

    gp = [-1/sqrt(3),  1/sqrt(3)];
    Ke = zeros(8,8);
    Me = zeros(8,8);

    xiN  = [-1  1  1 -1];
    etaN = [-1 -1  1  1];

    detJ = (dx*dy)/4;
    invJ = [2/dx 0; 0 2/dy];

    for i = 1:2
        for j = 1:2
            xi  = gp(i);
            eta = gp(j);

            N = 0.25 * (1 + xi*xiN) .* (1 + eta*etaN);
            dN_dxi  = 0.25 * xiN  .* (1 + eta*etaN);
            dN_deta = 0.25 * etaN .* (1 + xi*xiN);

            grads = invJ * [dN_dxi; dN_deta];
            dN_dx = grads(1,:);
            dN_dy = grads(2,:);

            B = zeros(3,8);
            for a = 1:4
                B(1,2*a-1) = dN_dx(a);
                B(2,2*a)   = dN_dy(a);
                B(3,2*a-1) = dN_dy(a);
                B(3,2*a)   = dN_dx(a);
            end

            Nmat = zeros(2,8);
            for a = 1:4
                Nmat(1,2*a-1) = N(a);
                Nmat(2,2*a)   = N(a);
            end

            Ke = Ke + (B' * C * B) * (t * detJ);
            Me = Me + (Nmat' * Nmat) * (rho * t * detJ);
        end
    end

    % 2x2 Gauss weights are 1*1, already implied above.
end

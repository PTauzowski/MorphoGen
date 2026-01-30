function [omega_best, xPhys_best] = topFreqOptimization_MMA_bound(L,H,nelx,nely,volfrac,penal,rmin,maxiter,supportType,J)

% ==============================================================
% Frequency maximization via MMA (BOUND formulation)
% Maximize min_{j=1..J} lambda_j  using variable Eb = E/lambda_ref
% ==============================================================
%
% Robust version:
%  - shift-invert eigensolver
%  - scaled bound variable
%  - MMA shape hardening (NO asymptote crashes)
%
% Output:
%   omega_best  [rad/s]
%   xPhys_best  physical density field
%
% ==============================================================

if nargin < 10 || isempty(J), J = 2; end   % J=2 is safest

%% --- material -------------------------------------------------
E0   = 1e7;
Emin = 1e-3;
rho0 = 1.0;
rho_min = 1e-6;
nu   = 0.3;
t    = 1.0;

%% --- mesh -----------------------------------------------------
dx = L/nelx; dy = H/nely;
nEl  = nelx*nely;

nodeNrs = reshape(1:(nelx+1)*(nely+1), nely+1, nelx+1);
nDof    = 2*(nelx+1)*(nely+1);

fixed = buildSupports(supportType,nodeNrs);
free  = setdiff(1:nDof,fixed);

%% --- element matrices ----------------------------------------
[Ke0, M0] = q4_rect_KeM_planeStress(E0,nu,rho0,t,dx,dy);
Ke = Ke0/E0;
Me = M0/rho0;

%% --- filter ---------------------------------------------------
[dyf,dxf] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1);
h  = max(0,rmin-sqrt(dxf.^2+dyf.^2));
Hs = imfilter(ones(nely,nelx),h,'symmetric');

fwd  = @(x) imfilter(x,h,'symmetric')./Hs;
bwd  = @(g) imfilter(g./Hs,h,'symmetric');

%% --- projection ----------------------------------------------
eta = 0.5;
betaMax = 16;

%% --- assembly indices ----------------------------------------
cVec = 2*nodeNrs(1:nely,1:nelx)+1;
cVec = cVec(:);
cMat = [cVec cVec+1 cVec+2*nely+2 cVec+2*nely+3 ...
        cVec+2*nely cVec+2*nely+1 cVec-2 cVec-1];

[Il,Jl] = find(tril(ones(8)));
iK = reshape(cMat(:,Il)',[],1);
jK = reshape(cMat(:,Jl)',[],1);
Ke_l = Ke(sub2ind([8,8],Il,Jl));
Me_l = Me(sub2ind([8,8],Il,Jl));

%% --- MMA setup -----------------------------------------------
lambda_ref = 2e4;         % scaling reference
n = nEl + 1;
m = J + 1;

xmin = [zeros(nEl,1); 0];
xmax = [ones(nEl,1);  5];      % Eb <= 5

xval  = [volfrac*ones(nEl,1); 1];
xold1 = xval;
xold2 = xval;
low   = xmin;
upp   = xmax;

a0 = 1;
a  = zeros(m,1);
c  = 10*ones(m,1);
d  = 1e-3*ones(m,1);

omega_best = -inf;
xPhys_best = xval(1:nEl);

%% =================== OPT LOOP ================================
for it = 1:maxiter

    % --- MMA shape hardening (CRITICAL)
    xval=xval(:); xold1=xold1(:); xold2=xold2(:);
    xmin=xmin(:); xmax=xmax(:);
    low=low(:); upp=upp(:);

    beta = min(betaMax,2^floor((it-1)/40));

    x = xval(1:nEl);
    Eb = xval(end);

    % --- filter + projection
    xT = fwd(reshape(x,nely,nelx));
    [xPhysMat,dH] = heavisideProjection(xT,beta,eta);
    xPhys = xPhysMat(:);

    % --- assemble K,M
    Ee = Emin + (xPhys.^penal)*(E0-Emin);
    re = rho_min + xPhys*(rho0-rho_min);

    K = sparse(iK,jK,kron(Ee,Ke_l),nDof,nDof);
    M = sparse(iK,jK,kron(re,Me_l),nDof,nDof);
    K = K+K'-diag(diag(K));
    M = M+M'-diag(diag(M));

    Kf = K(free,free);
    Mf = M(free,free);

    % --- eigensolve (SHIFT-INVERT!)
    sigma = Eb*lambda_ref;
    opts.tol=1e-8; opts.maxit=300; opts.disp=0;
    [V,D] = eigs(Kf,Mf,J,sigma,opts);
    lam = sort(real(diag(D)));

    omega_cur = sqrt(lam(1));
    if omega_cur > omega_best
        omega_best = omega_cur;
        xPhys_best = xPhys;
    end

    % --- sensitivities
    dlam = zeros(nEl,J);
    for j=1:J
        v = V(:,j); v=v/sqrt(v'*(Mf*v));
        phi=zeros(nDof,1); phi(free)=v;
        pe=phi(cMat);
        dlam(:,j)= ...
            penal*(E0-Emin)*(xPhys.^(penal-1)).*sum((pe*Ke).*pe,2) ...
          - lam(j)*(rho0-rho_min).*sum((pe*Me).*pe,2);
    end

    % --- MMA objective
    f0 = -Eb;
    df0 = zeros(n,1); df0(end)=-1;

    % --- constraints
    fval=zeros(m,1); dfdx=zeros(m,n);

    for j=1:J
        fval(j)=Eb-lam(j)/lambda_ref;
        g = reshape(-dlam(:,j)/lambda_ref,nely,nelx).*dH;
        dfdx(j,1:nEl)=bwd(g)(:);
        dfdx(j,end)=1;
    end

    fval(end)=mean(xPhys)-volfrac;
    gv=reshape((1/nEl)*ones(nEl,1),nely,nelx).*dH;
    dfdx(end,1:nEl)=bwd(gv)(:);

    % --- MMA step
    [xnew,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,it,xval,xmin,xmax,xold1,xold2,...
               f0,df0,fval,dfdx,low,upp,a0,a,c,d);

    xold2=xold1; xold1=xval; xval=xnew;

    fprintf('It:%3d  beta:%2d  omega:%.3f  vol:%.3f  maxg:%+.2e\n',...
        it,beta,omega_cur,mean(xPhys),max(fval))

    imagesc(1-reshape(xPhys,nely,nelx)); axis equal off; drawnow
end

fprintf('\nBest omega = %.4f rad/s\n',omega_best)
end


% ---------------- helpers ---------------------------------------
function [xPhys, dH] = heavisideProjection(xTilde, beta, eta)
    denom = tanh(beta*eta) + tanh(beta*(1-eta));
    xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / denom;
    dH    = (beta * (1 - tanh(beta*(xTilde-eta)).^2)) / denom;
end

function fixedDOFs = buildSupports(supportType, nodeNrs)
    nely = size(nodeNrs,1)-1;
    leftNodes  = nodeNrs(:,1);
    rightNodes = nodeNrs(:,end);

    leftBot  = nodeNrs(1,1);
    rightBot = nodeNrs(1,end);

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
            fixedDOFs = [u(leftBot); v(leftBot); v(rightBot)];
        otherwise
            error("Unknown supportType '%s'. Use CF, CC, CH, HH.", s);
    end
    fixedDOFs = unique(fixedDOFs(:));
end

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
end

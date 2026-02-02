function [KE, ME] = fe_hex8(E, nu, rho, EL, EW, EH)
%FE_HEX8 Compute 8-node hexahedral element stiffness and mass (24x24).
% Uses standard 2x2x2 Gauss integration for an orthogonal brick of size
% EL (x), EW (y), EH (z). Matches Appendix G intent without relying on OCR
% coefficient tables.

    % Material matrix for 3D elasticity
    lam = E*nu/((1+nu)*(1-2*nu));
    mu  = E/(2*(1+nu));
    C = [lam+2*mu, lam     , lam     , 0 , 0 , 0;
         lam     , lam+2*mu, lam     , 0 , 0 , 0;
         lam     , lam     , lam+2*mu, 0 , 0 , 0;
         0       , 0       , 0       , mu, 0 , 0;
         0       , 0       , 0       , 0 , mu, 0;
         0       , 0       , 0       , 0 , 0 , mu];

    % Gauss points and weights
    gp = 1/sqrt(3);
    gauss = [-gp, -gp, -gp;
              gp, -gp, -gp;
              gp,  gp, -gp;
             -gp,  gp, -gp;
             -gp, -gp,  gp;
              gp, -gp,  gp;
              gp,  gp,  gp;
             -gp,  gp,  gp];
    w = ones(8,1);  % all weights =1 for 2-point rule

    % Precompute Jacobian factors for orthogonal brick
    dxdxi = EL/2; dydet = EW/2; dzdze = EH/2;
    detJ  = dxdxi * dydet * dzdze;
    invJ = diag([2/EL, 2/EW, 2/EH]); % maps dN/dxi -> dN/dx

    KE = zeros(24,24);
    ME = zeros(24,24);

    for k = 1:8
        xi  = gauss(k,1); eta = gauss(k,2); zeta = gauss(k,3);
        [N, dN_dxi] = shape_hex8(xi, eta, zeta);
        dN_dx = (invJ * dN_dxi)'; % 8x3

        B = zeros(6,24);
        for a = 1:8
            id = (a-1)*3 + (1:3);
            dNx = dN_dx(a,1); dNy = dN_dx(a,2); dNz = dN_dx(a,3);
            B(:, id) = [ dNx,   0,    0;
                         0,     dNy,  0;
                         0,     0,    dNz;
                         dNy,   dNx,  0;
                         0,     dNz,  dNy;
                         dNz,   0,    dNx ];
        end

        KE = KE + (B' * C * B) * detJ * w(k);
        Nmat = zeros(3,24);
        for a = 1:8
            id = (a-1)*3 + (1:3);
            Nmat(:, id) = N(a) * eye(3);
        end
        ME = ME + rho * (Nmat' * Nmat) * detJ * w(k);
    end

    % Symmetrize (guard numerical noise)
    KE = (KE + KE')/2;
    ME = (ME + ME')/2;
end

function [N, dN_dxi] = shape_hex8(xi, eta, zeta)
% Shape functions and derivatives w.r.t local coords (xi,eta,zeta)
    N = 1/8 * [(1-xi)*(1-eta)*(1-zeta);
               (1+xi)*(1-eta)*(1-zeta);
               (1+xi)*(1+eta)*(1-zeta);
               (1-xi)*(1+eta)*(1-zeta);
               (1-xi)*(1-eta)*(1+zeta);
               (1+xi)*(1-eta)*(1+zeta);
               (1+xi)*(1+eta)*(1+zeta);
               (1-xi)*(1+eta)*(1+zeta)];

    dN_dxi = 1/8 * [ -(1-eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1-eta);
                     (1-eta)*(1-zeta),  -(1+xi)*(1-zeta), -(1+xi)*(1-eta);
                     (1+eta)*(1-zeta),   (1+xi)*(1-zeta), -(1+xi)*(1+eta);
                    -(1+eta)*(1-zeta),   (1-xi)*(1-zeta), -(1-xi)*(1+eta);
                    -(1-eta)*(1+zeta),  -(1-xi)*(1+zeta),  (1-xi)*(1-eta);
                     (1-eta)*(1+zeta),  -(1+xi)*(1+zeta),  (1+xi)*(1-eta);
                     (1+eta)*(1+zeta),   (1+xi)*(1+zeta),  (1+xi)*(1+eta);
                    -(1+eta)*(1+zeta),   (1-xi)*(1+zeta),  (1-xi)*(1+eta)]';
end

classdef Frame3DSectionProps < Frame3D
% Frame3DSectionProps  Frame3D variant with independent effective section
%   stiffnesses instead of geometry-derived ones.
%
%   Replaces the single-E surrogate in Frame3D with six independent channel
%   stiffnesses EA, EIy, EIz, GJ, GAy, GAz supplied directly by the caller.
%   All transformation and result-extraction methods are inherited from Frame3D
%   unchanged; only computeLocalStifnessMatrix is overridden.
%
%   Constructor:
%     obj = Frame3DSectionProps(elems, EA, EIy, EIz, GJ, GAy, GAz)
%
%   Inputs
%     elems  - [nElems x 2] connectivity (inherited Frame3D format)
%     EA     - axial stiffness [N]
%     EIy    - bending stiffness about local y [N·m²]
%     EIz    - bending stiffness about local z [N·m²]
%     GJ     - torsional stiffness [N·m²]
%     GAy    - shear stiffness in local y (Timoshenko, with kappa) [N]
%     GAz    - shear stiffness in local z (Timoshenko, with kappa) [N]

    properties
        EA_eff
        EIy_eff
        EIz_eff
        GJ_eff
        GAy_eff
        GAz_eff
    end

    methods
        function obj = Frame3DSectionProps(elems, EA, EIy, EIz, GJ, GAy, GAz)
            % Pass dummy geometry to Frame3D base constructor.
            % The base-class fields (E, A, G, Jy, Jz, Ks, kappa) are set
            % but never reached because computeLocalStifnessMatrix is
            % overridden below.
            E_dummy  = 1;
            nu_dummy = 0.3;
            R_dummy  = 1;
            r_dummy  = 0;
            obj = obj@Frame3D(elems, E_dummy, nu_dummy, R_dummy, r_dummy);

            obj.EA_eff  = EA;
            obj.EIy_eff = EIy;
            obj.EIz_eff = EIz;
            obj.GJ_eff  = GJ;
            obj.GAy_eff = GAy;
            obj.GAz_eff = GAz;
        end

        function K = computeLocalStifnessMatrix(obj, nodes, varargin)
        % Override: use effective stiffness scalars directly.
        % Timoshenko beam; Phi = 12*EI / (GA * l^2).
            nelems = size(obj.elems, 1);
            dim    = 12;
            K = zeros(dim, dim, nelems);

            EA    = obj.EA_eff;
            EJy   = obj.EIy_eff;
            EJz   = obj.EIz_eff;
            GKs   = obj.GJ_eff;
            kGAy  = obj.GAy_eff;
            kGAz  = obj.GAz_eff;

            for k = 1:nelems
                l  = norm(nodes(obj.elems(k,2),:) - nodes(obj.elems(k,1),:));
                l2 = l^2;
                l3 = l^3;

                Phy = 12*EJy / (kGAy * l2);
                Phz = 12*EJz / (kGAz * l2);
                cy  = 1 / (1 + Phy);
                cz  = 1 / (1 + Phz);

                Ke = zeros(12, 12);

                % axial (DOFs 1,7)
                Ke(1,1)   =  EA/l;   Ke(7,7)   =  EA/l;
                Ke(1,7)   = -EA/l;   Ke(7,1)   = -EA/l;

                % torsion (DOFs 4,10)
                Ke(4,4)   =  GKs/l;  Ke(10,10) =  GKs/l;
                Ke(4,10)  = -GKs/l;  Ke(10,4)  = -GKs/l;

                % xy-plane bending (EJz): DOFs 2,6 / 8,12
                Ke(2,2)   =  12*EJz*cz/l3;
                Ke(8,8)   =  12*EJz*cz/l3;
                Ke(6,6)   =  (4+Phz)*EJz*cz/l;
                Ke(12,12) =  (4+Phz)*EJz*cz/l;
                Ke(6,2)   =   6*EJz*cz/l2;  Ke(2,6)   =   6*EJz*cz/l2;
                Ke(8,2)   = -12*EJz*cz/l3;  Ke(2,8)   = -12*EJz*cz/l3;
                Ke(8,6)   =  -6*EJz*cz/l2;  Ke(6,8)   =  -6*EJz*cz/l2;
                Ke(12,2)  =   6*EJz*cz/l2;  Ke(2,12)  =   6*EJz*cz/l2;
                Ke(12,6)  =  (2-Phz)*EJz*cz/l; Ke(6,12) =  (2-Phz)*EJz*cz/l;
                Ke(12,8)  =  -6*EJz*cz/l2;  Ke(8,12)  =  -6*EJz*cz/l2;

                % xz-plane bending (EJy): DOFs 3,5 / 9,11
                Ke(3,3)   =  12*EJy*cy/l3;
                Ke(9,9)   =  12*EJy*cy/l3;
                Ke(5,5)   =  (4+Phy)*EJy*cy/l;
                Ke(11,11) =  (4+Phy)*EJy*cy/l;
                Ke(5,3)   =  -6*EJy*cy/l2;  Ke(3,5)   =  -6*EJy*cy/l2;
                Ke(9,3)   = -12*EJy*cy/l3;  Ke(3,9)   = -12*EJy*cy/l3;
                Ke(9,5)   =   6*EJy*cy/l2;  Ke(5,9)   =   6*EJy*cy/l2;
                Ke(11,3)  =  -6*EJy*cy/l2;  Ke(3,11)  =  -6*EJy*cy/l2;
                Ke(11,5)  =  (2-Phy)*EJy*cy/l; Ke(5,11) =  (2-Phy)*EJy*cy/l;
                Ke(11,9)  =   6*EJy*cy/l2;  Ke(9,11)  =   6*EJy*cy/l2;

                K(:,:,k) = Ke;
            end
        end
    end
end

classdef ReferenceModuleSolidModel < handle
    % Half-segment solid model for reference-module SIMP topology optimization.
    %
    % Geometry convention:
    %   Flat base at z=0   = LOADED face (section resultants applied here).
    %   Inclined top face  = FIXED face  (joint constraint).
    %
    % Load case struct fields (all in the solid coordinate frame):
    %   lc.N   – axial force  (+z = tension)
    %   lc.Ty  – shear in +y
    %   lc.Tz  – shear in +x
    %   lc.My  – bending moment about y-axis   (positive: tension at x>xc)
    %   lc.Mz  – bending moment about x-axis   (signed per frame-to-solid mapping)
    %   lc.Ms  – torsion about z-axis (CCW positive when viewed from +z)
    %
    % Resultant sign check used by validateLoadApplication:
    %   Fx_sum = lc.Tz,  Fy_sum = lc.Ty,  Fz_sum = -lc.N
    %   My_sum = lc.My   (= -sum(Pnodal_z .* (x-xc)))
    %   Mz_sum = lc.Mz   (=  sum(Pnodal_z .* (y-yc)))
    %   Ms_sum = lc.Ms   (=  sum(Pnodal_y .* (x-xc) - Pnodal_x .* (y-yc)))

    properties
        mesh
        fe
        analysis
        loadedFaceSelector
        fixedFaceSelector
        loaded_node_ids   % indices into mesh.nodes for the flat-base (z≈0) nodes

        % Nominal section properties (annular ring)
        A     % pi*(R^2 - r^2)
        Iy    % pi*(R^4 - r^4)/4
        Iz    % same as Iy by symmetry
        J     % pi*(R^4 - r^4)/2

        % Centroid of the loaded face (computed from actual node positions)
        xc
        yc
    end

    methods
        function obj = ReferenceModuleSolidModel(E, nu, r, R, h, alpha_deg, res_th)
            % Build the half-segment solid model.
            %   E, nu      – material parameters
            %   r, R       – inner/outer radii
            %   h          – segment length
            %   alpha_deg  – inclination of the top face (degrees)
            %   res_th     – radial resolution (circumferential/axial scaled)

            alpha_rad = alpha_deg * pi / 180;
            resCirc   = round(2*pi*R / (R-r) * res_th);
            resLen    = max(1, round(h / (R-r) * res_th));

            ShapeFn = ShapeFunctionH8;
            obj.mesh = Mesh();
            obj.mesh.addManipulatorHalfSegment3D(r, R, h, alpha_rad, ...
                                                 res_th, resCirc, resLen, ...
                                                 ShapeFn.localNodes);

            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();

            obj.fe = SolidElasticElem(ShapeFn, obj.mesh.elems);
            obj.fe.props.h = 1;
            obj.fe.setMaterial(material);

            obj.analysis = LinearElasticityWeighted(obj.fe, obj.mesh, false);

            % flat base = findDownwardFacingNodes (outward normal ≈ -z)
            % inclined top = findUpwardFacingNodes (outward normal has +z component)
            down_ids = obj.mesh.findDownwardFacingNodes();
            up_ids   = obj.mesh.findUpwardFacingNodes();

            dn_flag = false(size(obj.mesh.nodes, 1), 1);
            up_flag = false(size(obj.mesh.nodes, 1), 1);
            dn_flag(down_ids) = true;
            up_flag(up_ids)   = true;

            obj.loadedFaceSelector = Selector(dn_flag);
            obj.fixedFaceSelector  = Selector(up_flag);
            obj.loaded_node_ids    = down_ids(:);

            % Fix the inclined top face (joint)
            obj.analysis.fixNodes(obj.fixedFaceSelector, ["ux" "uy" "uz"]);

            % Nominal section properties
            obj.A  = pi * (R^2 - r^2);
            obj.Iy = pi * (R^4 - r^4) / 4;
            obj.Iz = pi * (R^4 - r^4) / 4;
            obj.J  = pi * (R^4 - r^4) / 2;

            % Centroid from actual flat-base node coordinates
            xn = obj.mesh.nodes(obj.loaded_node_ids, :);
            obj.xc = mean(xn(:, 1));
            obj.yc = mean(xn(:, 2));
        end

        % ------------------------------------------------------------------
        function applyLoadCase(obj, lc)
            % Clear any previous load and apply all six section resultants.
            obj.analysis.clearCurrentLoad();
            if isfield(lc, 'N')  && lc.N  ~= 0,  obj.applyAxialForce(lc.N);  end
            if isfield(lc, 'Ty') && lc.Ty ~= 0,  obj.applyShearY(lc.Ty);     end
            if isfield(lc, 'Tz') && lc.Tz ~= 0,  obj.applyShearX(lc.Tz);     end
            if isfield(lc, 'My') && lc.My ~= 0,  obj.applyBendingY(lc.My);   end
            if isfield(lc, 'Mz') && lc.Mz ~= 0,  obj.applyBendingZ(lc.Mz);   end
            if isfield(lc, 'Ms') && lc.Ms ~= 0,  obj.applyTorsion(lc.Ms);     end
        end

        % ------------------------------------------------------------------
        function result = validateLoadApplication(obj, lc, tol)
            % Sum nodal loads on the loaded face and compare to targets.
            % Returns struct with per-component errors and a boolean 'passed'.
            if nargin < 3, tol = 0.02; end

            ids = obj.loaded_node_ids;
            Pn  = obj.analysis.Pnodal(ids, :);
            xn  = obj.mesh.nodes(ids, :);
            dx  = xn(:, 1) - obj.xc;
            dy  = xn(:, 2) - obj.yc;

            Fx_sum  =  sum(Pn(:, 1));
            Fy_sum  =  sum(Pn(:, 2));
            Fz_sum  =  sum(Pn(:, 3));
            My_sum  = -sum(Pn(:, 3) .* dx);
            Mz_sum  =  sum(Pn(:, 3) .* dy);
            Ms_sum  =  sum(Pn(:, 2) .* dx - Pn(:, 1) .* dy);

            N_t  = lcfield(lc,'N',0);
            Ty_t = lcfield(lc,'Ty',0);
            Tz_t = lcfield(lc,'Tz',0);
            My_t = lcfield(lc,'My',0);
            Mz_t = lcfield(lc,'Mz',0);
            Ms_t = lcfield(lc,'Ms',0);

            labels  = ["Fx(=Tz)","Fy(=Ty)","Fz(=-N)","My","Mz","Ms"];
            targets = [Tz_t, Ty_t, -N_t, My_t, Mz_t, Ms_t];
            actual  = [Fx_sum, Fy_sum, Fz_sum, My_sum, Mz_sum, Ms_sum];

            result.Fx     = Fx_sum;
            result.Fy     = Fy_sum;
            result.Fz     = Fz_sum;
            result.My     = My_sum;
            result.Mz     = Mz_sum;
            result.Ms     = Ms_sum;
            result.passed = true;

            fprintf('\n%-12s  %12s  %12s  %8s\n','Component','Target','Actual','Error%');
            fprintf('%s\n', repmat('-',1,50));
            for k = 1:6
                t = targets(k);  a = actual(k);
                if abs(t) > 1e-14
                    err_pct = abs(a-t)/abs(t)*100;
                else
                    err_pct = abs(a)*100;
                end
                ok = (err_pct <= tol*100);
                if ~ok, result.passed = false; end
                flag = ''; if ~ok, flag = '  <-- FAIL'; end
                fprintf('%-12s  %12.6g  %12.6g  %7.2f%%%s\n', ...
                    labels(k), t, a, err_pct, flag);
            end
            if result.passed
                fprintf('PASS: all resultants within %.1f%%\n', tol*100);
            else
                fprintf('FAIL: some resultants outside %.1f%% tolerance\n', tol*100);
            end
        end
    end

    % ======================================================================
    methods (Access = private)

        function applyAxialForce(obj, N)
            % σ_z = N/A;  traction on -z face: tz = -N/A
            % Target: sum(Pn_z) = -N
            fn = @(x) [zeros(size(x,1),2), repmat(-N/obj.A, size(x,1),1)];
            obj.applyAndScale1DOF(fn, 3, -N);
        end

        function applyShearY(obj, Ty)
            % Uniform shear: ty = Ty/A.  Target: sum(Pn_y) = Ty
            fn = @(x) [zeros(size(x,1),1), repmat(Ty/obj.A, size(x,1),1), zeros(size(x,1),1)];
            obj.applyAndScale1DOF(fn, 2, Ty);
        end

        function applyShearX(obj, Tz)
            % Uniform shear: tx = Tz/A.  Target: sum(Pn_x) = Tz
            fn = @(x) [repmat(Tz/obj.A, size(x,1),1), zeros(size(x,1),2)];
            obj.applyAndScale1DOF(fn, 1, Tz);
        end

        function applyBendingY(obj, My)
            % σ_z = My*(x-xc)/Iy;  traction: tz = -My*(x-xc)/Iy
            % Integral gives My_actual = -sum(Pn_z*(x-xc)); scale to My.
            xc = obj.xc;  Iy = obj.Iy;
            fn = @(x) [zeros(size(x,1),2), -My*(x(:,1)-xc)/Iy];
            Pn_before = obj.analysis.Pnodal;
            obj.analysis.elementLoadSurfaceIntegral("global", ...
                obj.loadedFaceSelector, ["ux" "uy" "uz"], fn);
            ids = obj.loaded_node_ids;
            dPz = obj.analysis.Pnodal(ids,3) - Pn_before(ids,3);
            xn  = obj.mesh.nodes(ids,:);
            My_actual = -sum(dPz .* (xn(:,1) - xc));
            if abs(My_actual) > eps
                obj.analysis.Pnodal(ids,3) = Pn_before(ids,3) + dPz*(My/My_actual);
            end
        end

        function applyBendingZ(obj, Mz)
            % Convention: Mz causes σ_z = -Mz*(y-yc)/Iz → tz = Mz*(y-yc)/Iz
            % Integral gives Mz_actual = sum(Pn_z*(y-yc)); scale to Mz.
            yc = obj.yc;  Iz = obj.Iz;
            fn = @(x) [zeros(size(x,1),2), Mz*(x(:,2)-yc)/Iz];
            Pn_before = obj.analysis.Pnodal;
            obj.analysis.elementLoadSurfaceIntegral("global", ...
                obj.loadedFaceSelector, ["ux" "uy" "uz"], fn);
            ids = obj.loaded_node_ids;
            dPz = obj.analysis.Pnodal(ids,3) - Pn_before(ids,3);
            xn  = obj.mesh.nodes(ids,:);
            Mz_actual = sum(dPz .* (xn(:,2) - yc));
            if abs(Mz_actual) > eps
                obj.analysis.Pnodal(ids,3) = Pn_before(ids,3) + dPz*(Mz/Mz_actual);
            end
        end

        function applyTorsion(obj, Ms)
            % τ = Ms*r_loc/J in tangential direction.
            % tx = -Ms*(y-yc)/J,  ty = +Ms*(x-xc)/J
            % Ms_actual = sum(Pn_y*(x-xc) - Pn_x*(y-yc)); scale to Ms.
            xc = obj.xc;  yc = obj.yc;  J = obj.J;
            fn = @(x) Ms/J * [-(x(:,2)-yc), (x(:,1)-xc), zeros(size(x,1),1)];
            Pn_before = obj.analysis.Pnodal;
            obj.analysis.elementLoadSurfaceIntegral("global", ...
                obj.loadedFaceSelector, ["ux" "uy" "uz"], fn);
            ids = obj.loaded_node_ids;
            dPx = obj.analysis.Pnodal(ids,1) - Pn_before(ids,1);
            dPy = obj.analysis.Pnodal(ids,2) - Pn_before(ids,2);
            xn  = obj.mesh.nodes(ids,:);
            dx  = xn(:,1) - xc;
            dy  = xn(:,2) - yc;
            Ms_actual = sum(dPy.*dx - dPx.*dy);
            if abs(Ms_actual) > eps
                scale = Ms / Ms_actual;
                obj.analysis.Pnodal(ids,1) = Pn_before(ids,1) + dPx*scale;
                obj.analysis.Pnodal(ids,2) = Pn_before(ids,2) + dPy*scale;
            end
        end

        % ------------------------------------------------------------------
        function applyAndScale1DOF(obj, fn, dof, target)
            % Apply surface integral fn, then scale DOF dof on loaded face
            % so that sum(increment) == target.
            Pn_before = obj.analysis.Pnodal;
            obj.analysis.elementLoadSurfaceIntegral("global", ...
                obj.loadedFaceSelector, ["ux" "uy" "uz"], fn);
            ids  = obj.loaded_node_ids;
            dP   = obj.analysis.Pnodal(ids, dof) - Pn_before(ids, dof);
            fsum = sum(dP);
            if abs(fsum) > eps
                obj.analysis.Pnodal(ids, dof) = Pn_before(ids, dof) + dP*(target/fsum);
            end
        end
    end
end

% --------------------------------------------------------------------------
function v = lcfield(s, fname, dflt)
    if isfield(s, fname), v = s.(fname); else, v = dflt; end
end

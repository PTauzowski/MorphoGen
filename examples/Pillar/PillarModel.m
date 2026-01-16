classdef PillarModel < ModelLinear

    properties
        top_R
        pillar_layers, ground_layers
        pillar_res, ground_res
        pillar_chem, ground_chem
        int_th
        bank_h = 0;      % wysokość banku (<= suma "cap" warstw), ustaw w konstruktorze albo tu
        bank_r0 = 0;  
        sf
    end

    methods
        function obj = PillarModel( top_R, pillar_layers, ground_layers, pillar_res, ground_res, pillar_chem, ground_chem, int_th, sf )
            obj.top_R = top_R;
            obj.pillar_layers = pillar_layers;
            obj.ground_layers = ground_layers;
            obj.pillar_res = pillar_res;
            obj.ground_res = ground_res;
            obj.pillar_chem = pillar_chem;
            obj.ground_chem = ground_chem;
            obj.int_th = int_th;
            obj.sf = sf;
            obj.generateMesh();
        end

        function obj = generateMesh(obj)
            obj.mesh = Mesh();
            pillar_height = sum(obj.pillar_layers);
            ground_depth  = sum(obj.ground_layers);
            nrXY = 6;

            Rbase = obj.top_R + pillar_height * tan(5*pi/180);
            Rtop  = obj.top_R;

                        % ground matches pillar base
            obj.mesh.addLayeredQuarterCylinder([0, 0, -ground_depth], ...
                Rbase, obj.ground_layers, obj.ground_res, nrXY, true, obj.int_th, obj.sf.localNodes);
            
            % --- interface layer (between ground top and pillar bottom) ---
            % thickness = obj.int_th, placed as a thin layer just below z=0
            int_layers = obj.int_th;                       % single layer
            int_res    = max(1, obj.ground_res(end));      % or: max(1, obj.pillar_res(1)) / or dedicated parameter
            obj.mesh.addLayeredQuarterCylinder([0, 0, -obj.int_th], ...
                Rbase, int_layers, int_res, nrXY, true, obj.int_th, obj.sf.localNodes);
            
            % pillar then taper upwards
            obj.mesh.addLayeredQuarterCylinder([0, 0, 0], ...
                Rbase, obj.pillar_layers, obj.pillar_res, nrXY, true, obj.int_th, obj.sf.localNodes);
            
            obj.mesh.nodes = obj.mesh.coneTransformationX(pillar_height, Rbase, Rtop, obj.mesh.nodes);


            % ---- ring3D (pipe) under top ground layer ----
            z_bot = -obj.ground_layers(end);
            z_top = 0.0;
            Rin   = Rbase;
            Rout  = 1.2*Rbase;
            bank_hres = 5;


            obj.mesh.addPipe3D([0,0,z_bot], Rin, Rout, 0, 90, z_bot, z_top, bank_hres, 12, obj.ground_res(end), obj.sf.localNodes);

            % ---- depression applied to ring interior ONLY (boundary nodes fixed) ----
            % depth: how deep the depression is at the TOP surface (z=z_top), positive value.
            depth = 0.05 * abs(z_bot - z_top);   % <- tune (e.g. 0.2..0.6)*layer_thickness

            obj.mesh.nodes = PillarModel.applyRingDepressionZ( ...
                obj.mesh.nodes, [0 0], Rin, Rout, z_bot, z_top, depth );

            % ---- MISSING: deeper ground annulus behind trench (all layers except last) ----
            if numel(obj.ground_layers) >= 2
                deep_layers = obj.ground_layers(1:end-1);
                deep_res    = obj.ground_res(1:end-1);
            
                z0_deep = -ground_depth;          % bottom of full ground stack
                % top of deep part is exactly z_bot = -obj.ground_layers(end)
            
                obj.mesh = PillarModel.addLayeredPipe3D_noTopInterface( ...
                    obj.mesh, [0 0 0], Rin, Rout, 0, 90, ...
                    z0_deep, deep_layers, deep_res, true, obj.int_th, ...
                    bank_hres, 12, obj.sf.localNodes );
            end


            Rbank =  1.5 * Rout;

            obj.mesh = PillarModel.addLayeredPipe3D( ...
                obj.mesh, [0 0 0], Rout, Rbank, 0, 90, ...
                -ground_depth, obj.ground_layers, obj.ground_res, true, obj.int_th, ...
                bank_hres, 12, obj.sf.localNodes );


            obj_sf = ShapeFunctionH8();

            angle=50;


            x = [ Rout, 0,  z_top; Rbank, 0,  z_top; Rout + obj.pillar_layers(1) * tan(deg2rad(angle)), 0,  obj.pillar_layers(1) + z_top; Rbank,  0, obj.pillar_layers(1) + z_top; ...
                  Rout, deg2rad(90), z_top; Rbank, deg2rad(90), z_top; Rout + obj.pillar_layers(1) * tan(deg2rad(angle)), deg2rad(90), obj.pillar_layers(1) + z_top; Rbank, deg2rad(90), obj.pillar_layers(1) + z_top ];

            mesh1 = Mesh();
            mesh1.addRing3D( [0 0 0], obj_sf, x, bank_hres, obj.pillar_res(1),12 , obj.sf.localNodes);

            obj.mesh.mergeMesh(mesh1);

            Rin   = Rbase;
            Rout  = Rout + obj.pillar_layers(1) * tan(deg2rad(angle));

            z_top = obj.pillar_layers(1) + z_top;
            angle=60;
            layer2=2*obj.pillar_layers(1);

            x = [ Rout, 0,  z_top; Rbank, 0,  z_top; Rout + layer2 * tan(deg2rad(angle)), 0,  layer2 + z_top; Rbank,  0, layer2 + z_top; ...
                  Rout, deg2rad(90), z_top; Rbank, deg2rad(90), z_top; Rout + layer2 * tan(deg2rad(angle)), deg2rad(90), layer2 + z_top; Rbank, deg2rad(90), layer2 + z_top ];

            mesh1 = Mesh();
            mesh1.addRing3D( [0 0 0], obj_sf, x, bank_hres, 2, 12 , obj.sf.localNodes);

            obj.mesh.mergeMesh(mesh1);

            % ---- final circle -> square tile transition (preserve layering) ----
            Rtr   = Rbank;        % end of your cylindrical ground (already working)
            Rtile = 1.2*Rbank;    % size of square tile (tune)
            nThetaSeg = 12;        % 6..16 is typical

            obj.mesh = PillarModel.addLayeredTransitionCircleToSquare( ...
                obj.mesh, Rtr, Rtile, -ground_depth, [ obj.ground_layers obj.pillar_layers(1) layer2 ], [obj.ground_res obj.pillar_res(1) 2 ], ...
                true, obj.int_th, nThetaSeg, bank_hres, obj.sf );

           
         

        end
    end

    methods(Static)

        

        function nodes = warpCircleToSquareTileQuarter(nodes, x0_xy, Rtr, Rout)
            % Warp outer region from circular boundary to square tile (quarter domain),
            % preserving layered structure by keeping z unchanged.
            %
            % - No change for r <= Rtr
            % - Full square mapping at r >= Rout
            % - Smooth blend in between
        
            x = nodes(:,1) - x0_xy(1);
            y = nodes(:,2) - x0_xy(2);
        
            r = hypot(x,y);
            epsR = 1e-12 * max(1, Rout);
        
            mask = (r >= Rtr - epsR) & (r <= Rout + epsR) & (r > epsR);
            if ~any(mask), return; end
        
            idx = find(mask);
        
            rx = x(idx);
            ry = y(idx);
            rr = r(idx);
        
            ux = rx ./ rr;
            uy = ry ./ rr;
        
            % square intersection scale along direction u
            s = 1 ./ max(ux, uy);           % quarter domain: ux,uy >= 0
        
            % mapped coordinates on "square level set"
            x_sq = rr .* s .* ux;
            y_sq = rr .* s .* uy;
        
            % blend weight based on radius
            t = (rr - Rtr) ./ (Rout - Rtr);
            t = max(0, min(1, t));
        
            % smoothstep
            w = t .* t .* (3 - 2 .* t);
        
            % apply blend
            x_new = (1 - w) .* rx + w .* x_sq;
            y_new = (1 - w) .* ry + w .* y_sq;
        
            % freeze exact boundaries (avoid moving seam nodes numerically)
            isRtr  = abs(rr - Rtr)  <= 100*epsR;
            isRout = abs(rr - Rout) <= 100*epsR;
            x_new(isRtr | isRout) = rx(isRtr | isRout);
            y_new(isRtr | isRout) = ry(isRtr | isRout);
        
            nodes(idx,1) = x_new + x0_xy(1);
            nodes(idx,2) = y_new + x0_xy(2);
        end


        function mesh = addLayeredPipe3D_noTopInterface(mesh, x0, r0, r1, al1, al2, z0, h, nrZ, add_interface, i_th, nr, nc, lnodes)
            % Build layered pipe from z0 upwards, with optional interfaces BETWEEN layers,
            % but WITHOUT adding an interface slab above the last layer.
        
            h   = h(:);
            nrZ = nrZ(:);
        
            L = numel(h);
            if numel(nrZ) ~= L
                error('addLayeredPipe3D_noTopInterface: h and nrZ must have same length.');
            end
        
            add_interface = logical(add_interface);
        
            % effective heights like addLayeredQuarterCylinder (only internal interfaces)
            if add_interface
                if ~(isscalar(i_th) && i_th > 0)
                    error('addLayeredPipe3D_noTopInterface: i_th must be > 0.');
                end
                h_eff = h;
                for k = 1:L
                    if k > 1, h_eff(k) = h_eff(k) - i_th/2; end
                    if k < L, h_eff(k) = h_eff(k) - i_th/2; end
                end
                if any(h_eff <= 0)
                    error('addLayeredPipe3D_noTopInterface: i_th too large -> some h_eff <= 0.');
                end
            else
                h_eff = h;
            end
        
            zoff = 0.0;
        
            for k = 1:L
                % main layer
                z1 = z0 + zoff;
                z2 = z1 + h_eff(k);
                mesh.addPipe3D(x0, r0, r1, al1, al2, z1, z2, nr, nc, nrZ(k), lnodes);
                zoff = zoff + h_eff(k);
        
                % interface between k and k+1 only (no interface after the last one)
                if add_interface && (k < L)
                    z1 = z0 + zoff;
                    z2 = z1 + i_th;
                    mesh.addPipe3D(x0, r0, r1, al1, al2, z1, z2, nr, nc, 1, lnodes);
                    zoff = zoff + i_th;
                end
            end
        end


        function mesh = addLayeredPipe3D(mesh, x0, r0, r1, al1, al2, z0, h, nrZ, add_interface, i_th, nr, nc, lnodes)
            % Build annular quarter-pipe with the SAME layered structure as addLayeredQuarterCylinder.
            % z0 is the bottom (e.g. -ground_depth). Total height = sum(h).
        
            h   = h(:);
            nrZ = nrZ(:);

            L = numel(h);
            if numel(nrZ) ~= L
                error('addLayeredPipe3D: h and nrZ must have the same length.');
            end
        
            add_interface = logical(add_interface);
        
            if add_interface
                if ~(isscalar(i_th) && i_th > 0)
                    error('addLayeredPipe3D: i_th must be > 0 when add_interface==true.');
                end
                if L < 2
                    add_interface = false;
                end
            end
        
            % effective heights like Mesh.addLayeredQuarterCylinder
            if add_interface
                h_eff = h;
                for k = 1:L
                    if k > 1, h_eff(k) = h_eff(k) - i_th/2; end
                    if k < L, h_eff(k) = h_eff(k) - i_th/2; end
                end
                if any(h_eff <= 0)
                    error('addLayeredPipe3D: i_th too large -> some h_eff <= 0.');
                end
            else
                h_eff = h;
            end
        
            zoff = 0.0;
        
            for k = 1:L
                % main layer k
                z1 = z0 + zoff;
                z2 = z1 + h_eff(k);
                mesh.addPipe3D(x0, r0, r1, al1, al2, z1, z2, nr, nc, nrZ(k), lnodes);
                zoff = zoff + h_eff(k);
        
                % interface between k and k+1 (nz=1 like in your quarter-cylinder)
                if add_interface && (k < L)
                    z1 = z0 + zoff;
                    z2 = z1 + i_th;
                    mesh.addPipe3D(x0, r0, r1, al1, al2, z1, z2, nr, nc, 1, lnodes);
                    zoff = zoff + i_th;
                end
            end
        
            targetH = sum(h);
            if abs(zoff - targetH) > 1e-10 * max(1, targetH)
                error('addLayeredPipe3D: height mismatch built=%.15g expected=%.15g.', zoff, targetH);
            end
        end

        function nodes = applyRingDepressionZ(nodes, x0_xy, Rin, Rout, z_bot, z_top, depth)
            % Apply smooth "bowl" depression on a cylindrical ring volume.
            % Boundary is NOT moved:
            %   - r == Rin  (inner boundary)
            %   - r == Rout (outer boundary)
            %   - z == z_bot (bottom of the ring)
            %
            % Deformation acts strongest at z_top and fades to 0 at z_bot.

            x = nodes(:,1) - x0_xy(1);
            y = nodes(:,2) - x0_xy(2);
            z = nodes(:,3);

            r = hypot(x,y);

            % select ring volume (quarter is okay; we just use r,z)
            epsR = 1e-9 * max(1, Rout);
            epsZ = 1e-9 * max(1, abs(z_top - z_bot));

            inR = (r >= Rin - epsR) & (r <= Rout + epsR);
            inZ = (z >= z_bot - epsZ) & (z <= z_top + epsZ);
            mask = inR & inZ;

            if ~any(mask)
                return;
            end

            % normalized radius in [0,1]
            t = (r(mask) - Rin) ./ (Rout - Rin);
            t = max(0, min(1, t));

            % radial bowl shape: 0 at t=0 and t=1, min at t=0.5
            bowl = 4 .* t .* (1 - t);     % in [0,1]

            % normalized height in [0,1] (0 at bottom, 1 at top)
            s = (z(mask) - z_bot) ./ (z_top - z_bot);
            s = max(0, min(1, s));

            % smoothstep to avoid kinks
            s = s.*s.*(3 - 2.*s);

            dz = -depth .* bowl .* s;

            % Freeze boundary nodes explicitly:
            isInner = abs(r(mask) - Rin)  <= 10*epsR;
            isOuter = abs(r(mask) - Rout) <= 10*epsR;
            isBot   = abs(z(mask) - z_bot) <= 10*epsZ;

            dz(isInner | isOuter | isBot) = 0;

            nodes(mask,3) = nodes(mask,3) + dz;
        end

        function mesh = addLayeredCapBank(mesh, r0, r1, th0, th1, z0, bank_h, layer_th, layer_res, nr, nt, localNodes)
            % Tworzy "cap bank" jako wiele cienkich pipe'ów (po jednej warstwie),
            % w zakresie r=[r0..r1], z=[0..bank_h] (nad gruntem).
            %
            % layer_th i layer_res takie jak pillar (warstwowanie identyczne).
    
            z = z0;
            for k = 1:numel(layer_th)
                if z >= z0 + bank_h
                    break;
                end
                dz = layer_th(k);
                if z + dz > z0 + bank_h
                    dz = (z0 + bank_h) - z;   % utnij ostatnią warstwę
                end
                nz = layer_res(min(k, numel(layer_res)));
                if isempty(nz) || nz < 1
                    nz = 1;
                end
    
                mesh.addPipe3D([0,0,z], r0, r1, th0, th1, z, z+dz, nr, nt, nz, localNodes);
                z = z + dz;
            end
        end
    
        function nodes = warpBankRamp(nodes, r0, r1, z_bot, z_top)
            x = nodes(:,1);
            y = nodes(:,2);
            z = nodes(:,3);
            r = hypot(x,y);
        
            epsR = 1e-9 * max(1, r1);
            epsZ = 1e-9 * max(1, abs(z_top - z_bot));
        
            inR = (r >= r0 - epsR) & (r <= r1 + epsR);
            inZ = (z >= z_bot - epsZ) & (z <= z_top + epsZ);
            mask = inR & inZ;
        
            if ~any(mask)
                return;
            end
        
            idx = find(mask);          % global indices of affected nodes
        
            t = (r(idx) - r0) ./ (r1 - r0);
            t = max(0, min(1, t));
        
            % smoothstep in [0,1]
            phi = t.*t.*(3 - 2.*t);
        
            s = (z(idx) - z_bot) ./ (z_top - z_bot);
            s = max(0, min(1, s));
            s = s.^1.6;  % tune
        
            z_target = z_top .* phi;
            z_new = z(idx) + s .* (z_target - z(idx));
        
            % Freeze boundaries (indices are LOCAL to idx)
            isR0  = abs(r(idx) - r0) <= 10*epsR;
            isR1  = abs(r(idx) - r1) <= 10*epsR;
            isBot = abs(z(idx) - z_bot) <= 10*epsZ;
            freeze = isR0 | isR1 | isBot;
        
            % restore original z for frozen nodes
            z_new(freeze) = z(idx(freeze));
        
            nodes(idx,3) = z_new;
        end

        function mesh = addLayeredTransitionCircleToSquare(mesh, Rin, Rtile, z0, h, nrZ, add_interface, i_th, nThetaSeg, nr, sf)
            % Builds a circle->square transition band for the GROUND:
            %   inner boundary: circle of radius Rin (conforming with existing cyl)
            %   outer boundary: square tile with half-width Rtile
            % Layering in Z follows h, nrZ, with optional interface slabs (like your ground).
        
            h   = h(:);
            nrZ = nrZ(:);
            L = numel(h);
            if numel(nrZ) ~= L
                error('Transition: h and nrZ must match.');
            end
        
            % effective heights like your layered generators
            if add_interface && L >= 2
                h_eff = h;
                for k = 1:L
                    if k > 1, h_eff(k) = h_eff(k) - i_th/2; end
                    if k < L, h_eff(k) = h_eff(k) - i_th/2; end
                end
            else
                h_eff = h;
                add_interface = false;
            end
        
            zoff = 0.0;
        
            th_edges = linspace(0, pi/2, nThetaSeg+1);
        
            for k = 1:L
                % main layer k
                z1 = z0 + zoff;
                z2 = z1 + h_eff(k);
        
                mesh = PillarModel.addTransitionThetaSectors(mesh, Rin, Rtile, th_edges, z1, z2, nr, nrZ(k), sf);
        
                zoff = zoff + h_eff(k);
        
                % interface slab between layers
                if add_interface && (k < L)
                    z1 = z0 + zoff;
                    z2 = z1 + i_th;
        
                    % interface uses nz=1 (consistent with your other generators)
                    mesh = PillarModel.addTransitionThetaSectors(mesh, Rin, Rtile, th_edges, z1, z2, nr, 1, sf);
        
                    zoff = zoff + i_th;
                end
            end
        end
        
        function mesh = addTransitionThetaSectors(mesh, Rin, Rtile, th_edges, z1, z2, nr, nz, sf)
            % Build multiple addRing3D blocks, one per theta-sector.
            % Inner radius is constant (Rin). Outer radius varies with theta to form square tile.
        
            obj_sf = ShapeFunctionH8();  % if your addRing3D needs it; MUST be compatible with sf.localNodes in merge
        
            for i = 1:(numel(th_edges)-1)
                thA = th_edges(i);
                thB = th_edges(i+1);
        
                RoutA = Rtile / max(cos(thA), sin(thA));
                RoutB = Rtile / max(cos(thB), sin(thB));
        
                % 8 corners in cylindrical coords [R, theta, z]
                % (inner/outer) x (theta A/B) x (z1/z2)
                x = [
                    Rin,  thA, z1;
                    RoutA,thA, z1;
                    Rin,  thA, z2;
                    RoutA,thA, z2;
        
                    Rin,  thB, z1;
                    RoutB,thB, z1;
                    Rin,  thB, z2;
                    RoutB,thB, z2;
                ];
        
                m = Mesh();
                % NOTE: adjust argument order to match YOUR Mesh.addRing3D signature.
                % Your earlier call suggests: addRing3D(x0, obj_sf, x, nr, nz, nc, localNodes)
                m.addRing3D([0 0 0], obj_sf, x, nr, nz, 1, sf.localNodes);
        
                mesh.mergeMesh(m);
            end
        end


    end
end

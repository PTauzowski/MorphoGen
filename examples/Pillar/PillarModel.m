classdef PillarModel < ModelLinear

    properties
        top_R
        pillar_layers, ground_layers, ground_layers_eff
        pillar_res, ground_res
        pillar_chem, ground_chem
        int_th
        bank_h = 0;      % wysokość banku (<= suma "cap" warstw), ustaw w konstruktorze albo tu
        bank_r0 = 0;  
        sf
        z_corners_coords
        z_coords
        z_tolerance
        z_offset
    end

    methods
        function obj = PillarModel( top_R, pillar_layers, ground_layers, pillar_res, ground_res, pillar_chem, ground_chem, int_th, sf, z_offset )
            obj.z_tolerance=1e-9;
            obj.top_R = top_R;
            obj.pillar_layers = pillar_layers;
            obj.ground_layers = ground_layers;
            obj.pillar_res = pillar_res;
            obj.ground_res = ground_res;
            obj.pillar_chem = pillar_chem;
            obj.ground_chem = ground_chem;
            obj.int_th = int_th;
            obj.sf = sf;
            obj.z_offset = z_offset;
            obj.generateMesh();
            obj.mesh.nodes(:,3) = obj.mesh.nodes(:,3) + obj.z_offset;
            obj.z_corners_coords = obj.z_corners_coords + obj.z_offset;
            obj.z_coords = obj.z_coords + obj.z_offset;
            obj.computeZNodalCoords();
        end


        function obj = generateMesh(obj)
            obj.mesh = Mesh();
            mergeTol = max(1e-6, obj.z_tolerance);
            obj.mesh.tolerance = mergeTol;
            pillar_height = sum(obj.pillar_layers);
            ground_depth  = sum(obj.ground_layers);
            nrXY = 4;

            Rbase = obj.top_R + pillar_height * tan(5*pi/180);
            Rtop  = obj.top_R;

            % ground matches pillar base, ends at z = -int_th/2 (bottom of centered interface)
            obj.ground_layers_eff = obj.ground_layers;
            if obj.int_th > 0
                obj.ground_layers_eff(end) = obj.ground_layers_eff(end) - obj.int_th/2;
                if obj.ground_layers_eff(end) <= 0
                    error('int_th too large: top ground layer would become <= 0.');
                end
            end
            
            obj.mesh.addLayeredQuarterCylinder([0, 0, -ground_depth], ...
                Rbase, obj.ground_layers_eff, obj.ground_res, nrXY, true, obj.int_th, obj.sf.localNodes);

            
            % --- interface layer (between ground top and pillar bottom) ---
            % Interface is CENTERED at z=0, spanning z in [-int_th/2, +int_th/2]
            % This takes int_th/2 from ground and int_th/2 from pillar nominal regions
            int_layers = obj.int_th;                       % single layer
            int_res    = 1;                                % interface layers have exactly 1 element in z-direction
            obj.mesh.addLayeredQuarterCylinder([0, 0, -obj.int_th/2], ...
                Rbase, int_layers, int_res, nrXY, true, obj.int_th, obj.sf.localNodes);

            % pillar starts at top of interface (z = +int_th/2)
            pillar_layers_eff = obj.pillar_layers;
            if obj.int_th > 0
                pillar_layers_eff(1) = pillar_layers_eff(1) - obj.int_th/2;
            end
            
            obj.mesh.addLayeredQuarterCylinder([0, 0, obj.int_th/2], ...
                Rbase, pillar_layers_eff, obj.pillar_res, nrXY, true, obj.int_th, obj.sf.localNodes);

            
            obj.mesh.nodes = obj.mesh.coneTransformationX(pillar_height, Rbase, Rtop, obj.mesh.nodes);


            % ---- ring3D (pipe) under top ground layer ----
            z0_ground = -ground_depth; % aligned with core ground stack
            z_top = -obj.int_th/2;                    % top of ground (below top interface)
            Rin   = Rbase;
            Rout  = 1.2*Rbase;
            bank_hres = 5;
  
            % Build full layered annulus with the same layering logic as the core.
            obj.mesh = PillarModel.addLayeredPipe3D( ...
                obj.mesh, [0 0 0], Rin, Rout, 0, 90, ...
                z0_ground, obj.ground_layers_eff, obj.ground_res, true, obj.int_th, ...
                bank_hres, nrXY*2, obj.sf.localNodes );

            % ---- depression applied to ring interior ONLY (boundary nodes fixed) ----
            % depth: how deep the depression is at the TOP surface (z=z_top), positive value.
            % Cap depth so top nodes do not cross (or nearly collide with) the
            % first interior z-level in the upper ground layer.
            depthReq = 0.05 * obj.ground_layers_eff(end);   % <- tune base ratio

            sfStep = max(1, numel(unique(obj.sf.localNodes(:,3))) - 1); % H8->1, H27->2
            topLayerEff = obj.ground_layers_eff(end);
            if obj.int_th > 0 && numel(obj.ground_layers_eff) >= 2
                topLayerEff = topLayerEff - obj.int_th/2;
            end
            dzTopNode = topLayerEff / max(1, sfStep * obj.ground_res(end));
            minClearance = max(obj.int_th, 1e-3);
            depthMax = max(0, dzTopNode - minClearance);
            depth = min(depthReq, depthMax);

            if depth < depthReq
                warning('PillarModel:DepressionDepthCapped', ...
                    ['Requested ring depression depth %.6g exceeds safe limit %.6g ' ...
                     '(top nodal dz %.6g, clearance %.6g). Using capped value.'], ...
                    depthReq, depthMax, dzTopNode, minClearance);
            end

            obj.mesh.nodes = PillarModel.applyRingDepressionZ( ...
                obj.mesh.nodes, [0 0], Rin, Rout, z0_ground, z_top, depth );

            


            Rbank =  1.5 * Rout;

            obj.mesh = PillarModel.addLayeredPipe3D( ...
                obj.mesh, [0 0 0], Rout, Rbank, 0, 90, ...
                z0_ground, obj.ground_layers_eff, obj.ground_res, true, obj.int_th, ...
                bank_hres, nrXY*2, obj.sf.localNodes );


            obj_sf = ShapeFunctionH8();

            angle=50;


            x = [ Rout, 0,  z_top;  Rbank, 0,  z_top;  Rout, deg2rad(90), z_top;  Rbank, deg2rad(90), z_top; ...
Rout + obj.pillar_layers(1) * tan(deg2rad(angle)), 0,  obj.pillar_layers(1) + z_top;  Rbank,  0, obj.pillar_layers(1) + z_top; Rout + obj.pillar_layers(1) * tan(deg2rad(angle)), deg2rad(90), obj.pillar_layers(1) + z_top; 
Rbank, deg2rad(90), obj.pillar_layers(1) + z_top  ];

            mesh1 = Mesh();
            mesh1.tolerance = mergeTol;
            mesh1.addRing3D( [0 0 0], obj_sf, x, bank_hres, nrXY*2, obj.pillar_res(1) , obj.sf.localNodes);

            obj.mesh.mergeMesh(mesh1);

            Rin   = Rbase;
            Rout  = Rout + obj.pillar_layers(1) * tan(deg2rad(angle));

            z_top = obj.pillar_layers(1) + z_top;
            angle=60;
            layer2=2*obj.pillar_layers(1);

            x = [ Rout, 0,  z_top;  Rbank, 0,  z_top; Rout, deg2rad(90), z_top; Rbank, deg2rad(90), z_top; ...
Rout + layer2 * tan(deg2rad(angle)), 0,  layer2 + z_top;  Rbank,  0, layer2 + z_top; Rout + layer2 * tan(deg2rad(angle)), deg2rad(90), layer2 + z_top; Rbank, deg2rad(90), layer2 + z_top ];

            mesh1 = Mesh();
            mesh1.tolerance = mergeTol;
            mesh1.addRing3D( [0 0 0], obj_sf, x, bank_hres,  nrXY*2, 2 , obj.sf.localNodes);

            obj.mesh.mergeMesh(mesh1);

            % ---- final circle -> square tile transition (preserve layering) ----
            Rtr   = Rbank;        % end of your cylindrical ground (already working)
            Rtile = 1.2*Rbank;    % size of square tile (tune)
            nThetaSeg = nrXY*2;        % 6..16 is typical

            % Ground transition: keep EXACT same ground layering as annulus/core.
            % This avoids z-level drift on the Rbank seam.
            obj.mesh = PillarModel.addLayeredTransitionCircleToSquare( ...
                obj.mesh, Rtr, Rtile, z0_ground, obj.ground_layers_eff, obj.ground_res, ...
                true, obj.int_th, nThetaSeg, bank_hres, obj.sf );

            % Upper bank transition (above ground top): keep horizontal
            % internal interfaces and preserve the bank-side geometry logic.
            z0_bank_top = -obj.int_th/2;
            obj.mesh = PillarModel.addLayeredTransitionCircleToSquare( ...
                obj.mesh, Rtr, Rtile, z0_bank_top, [obj.pillar_layers(1) layer2], [obj.pillar_res(1) 2], ...
                false, obj.int_th, nThetaSeg, bank_hres, obj.sf );
           % 
           % obj.smoothRectEdges(Rtile, Rtr, 12, 0.35);
         
           obj.mesh.nodes = obj.snapCriticalZCoordinates(obj.mesh.nodes);
           obj.snapAllCoordinates();
           weldStats = obj.mesh.weldNodes(mergeTol, true);
           fprintf('Global weld: nodes %d -> %d, elems %d -> %d (tol=%g)\n', ...
               weldStats.nodesBefore, weldStats.nodesAfter, ...
               weldStats.elemsBefore, weldStats.elemsAfter, mergeTol);

           % Try to correct residual orientation-only Jacobian issues.
           if size(obj.mesh.elems,2) == 27
               [badE0, ~] = obj.mesh.findNegativeJacobian(obj.sf, 1e-12);
               if ~isempty(badE0)
                   fixRes = obj.mesh.fixNegativeJacobianByRenumbering(obj.sf, 1e-12);
                   if fixRes.fixed_count > 0
                       fprintf('Jacobian renumbering: fixed %d elements, remaining bad: %d\n', ...
                           fixRes.fixed_count, numel(fixRes.badE_after));
                   end
               end
           end
           

        end

        
    
    function smoothRectEdges(obj, Rtile, Rtr, nIter, alpha)
        % smoothRectEdges  Straighten and smooth the OUTER square boundary.
        %
        % Rtile : target square half-width (outer boundary)
        % Rtr   : inner radius of transition band (where circle ends)
        % nIter : smoothing iterations (e.g. 5..20)  [optional, default 10]
        % alpha : relaxation (e.g. 0.2..0.6)        [optional, default 0.35]
        %
        % Effect:
        %   (1) Snap outer boundary nodes to x=Rtile or y=Rtile
        %   (2) Laplacian smooth XY in transition band only, keeping boundary fixed
        %       (z is never changed)

        if nargin < 4 || isempty(nIter),  nIter = 10; end
        if nargin < 5 || isempty(alpha),  alpha = 0.35; end

        nodes = obj.mesh.nodes;
        x = nodes(:,1); y = nodes(:,2); z = nodes(:,3);
        r = hypot(x,y);

        % -----------------------------
        % (A) SNAP outer boundary to exact square
        % -----------------------------
        tol = 1e-9 * max(1, Rtile);

        % outer boundary nodes: near max(x,y)=Rtile in the quarter
        outer = (x >= -tol) & (y >= -tol) & (max(x,y) >= (Rtile - 50*tol));

        % decide which side a node belongs to
        onX = outer & (x >= y);   % closer to x-wall
        onY = outer & (y >  x);   % closer to y-wall

        x(onX) = Rtile;
        y(onY) = Rtile;

        % make corner exact (optional but helps)
        corner = outer & (abs(x - Rtile) <= 100*tol) & (abs(y - Rtile) <= 100*tol);
        x(corner) = Rtile;
        y(corner) = Rtile;

        nodes(:,1) = x; nodes(:,2) = y;
        obj.mesh.nodes = nodes;

        % -----------------------------
        % (B) Build node adjacency from hex elements (for Laplacian smoothing)
        % -----------------------------
        elems = obj.mesh.elems;
        nN = size(nodes,1);
        adj = cell(nN,1);

        % connect nodes that share an element (cheap, good enough)
        for e = 1:size(elems,1)
            en = unique(elems(e,:));
            for ii = 1:numel(en)
                ni = en(ii);
                adj{ni} = [adj{ni}, en]; %#ok<AGROW>
            end
        end
        for i = 1:nN
            adj{i} = unique(adj{i}(adj{i} ~= i));
        end

        % -----------------------------
        % (C) Smooth only interior nodes of transition band (XY only)
        % -----------------------------
        nodes = obj.mesh.nodes;
        x = nodes(:,1); y = nodes(:,2); z = nodes(:,3); %#ok<NASGU>
        r = hypot(x,y);

        % band to smooth: between circle end and square boundary
        band = (r >= Rtr - 1e-9*Rtile) & (max(x,y) <= Rtile - 1e-9*Rtile);

        % fixed nodes: symmetry planes and snapped outer boundary
        fixed = (abs(x) <= tol) | (abs(y) <= tol) | outer;

        move = band & ~fixed;

        for it = 1:nIter
            xNew = x;
            yNew = y;

            idx = find(move);
            for k = 1:numel(idx)
                i = idx(k);
                nb = adj{i};
                if isempty(nb), continue; end
                mx = mean(x(nb));
                my = mean(y(nb));
                xNew(i) = (1-alpha)*x(i) + alpha*mx;
                yNew(i) = (1-alpha)*y(i) + alpha*my;
            end

            x = xNew;
            y = yNew;

            % keep fixed nodes exact
            x(onX) = Rtile;
            y(onY) = Rtile;
            x(abs(x) <= tol) = 0;
            y(abs(y) <= tol) = 0;
        end

        nodes(:,1) = x;
        nodes(:,2) = y;
        obj.mesh.nodes = nodes;
    end


    function [chem, lays_between_layers] = chemFromZ(obj, z)
        % chemFromZ  Map z-coordinates to chemistry with explicit interface slabs.
        %
        % Rules:
        % 1. Decide if in pillar (z>=0) or ground (z<0).
        % 2. Locate which segment (nominal vs interface).
        % 3. Return chemistry:
        %    - nominal segment → constant chem(layer)
        %    - interface slab → linear interpolation between chem_below and chem_above
        % 4. lays_between_layers:
        %    - true ONLY on planes that are boundaries between nominal and interface slabs
        %    - excluding global top and global bottom
        %    - false for interior points of interface slabs (including H27 mid-plane nodes)
        %
        % Note: The mesh generator builds ground layers from z=-depth upward, so
        % ground_layers(1) is at the bottom (most negative z) and ground_layers(end)
        % is at the top (closest to z=0). The special interface slab z in [-int_th, 0]
        % connects ground layer L (topmost) to pillar layer 1.

            % ---- validate ----
            if numel(obj.pillar_layers) ~= numel(obj.pillar_chem)
                error('pillar_layers (%d) and pillar_chem (%d) must match.', ...
                    numel(obj.pillar_layers), numel(obj.pillar_chem));
            end
            if numel(obj.ground_layers_eff) ~= numel(obj.ground_chem)
                error('ground_layers (%d) and ground_chem (%d) must match.', ...
                    numel(obj.ground_layers_eff), numel(obj.ground_chem));
            end

            z  = double(z);
            % Work in model-local z where the top interface is centered at 0.
            % FEAP export may shift all coordinates by obj.z_offset.
            z_local = z;
            if isprop(obj,'z_offset') && ~isempty(obj.z_offset)
                z_local = z - double(obj.z_offset);
            end
            sz = size(z);

            outIsCell = iscell(obj.pillar_chem) || iscell(obj.ground_chem);
            if outIsCell
                chem = cell(sz);
            else
                chem = zeros(sz);
            end
            lays_between_layers = false(sz);

            % strict tolerance for classification
            if isprop(obj,'z_tolerance') && ~isempty(obj.z_tolerance) && obj.z_tolerance > 0
                tol = obj.z_tolerance;
            else
                tol = 1e-9;
            end

            i_th = obj.int_th;

            % ---- helper: build 1D stack segments for nominal layers ----
            % Returns segment boundaries and interface planes in stack coordinate s
            % where s=0 is the start (bottom for pillar, TOP for ground after flip).
            function [segS0, segS1, segType, segK, segBelow, segAbove, planesBetween, totalH] = buildStack(h)
                h = h(:);
                L = numel(h);
                totalH = sum(h);

                if L == 0
                    segS0=[]; segS1=[]; segType=[]; segK=[]; segBelow=[]; segAbove=[]; planesBetween=[]; totalH=0;
                    return;
                end

                add_interface = (L >= 2) && (isscalar(i_th) && i_th > 0);

                % effective heights (must mirror mesh generator logic)
                if add_interface
                    h_eff = h;
                    for kk = 1:L
                        if kk > 1, h_eff(kk) = h_eff(kk) - i_th/2; end
                        if kk < L, h_eff(kk) = h_eff(kk) - i_th/2; end
                    end
                    if any(h_eff <= 0)
                        error('chemFromZ: interface thickness too large -> some effective heights <= 0.');
                    end
                else
                    h_eff = h;
                end

                segS0 = []; segS1 = []; segType = []; segK = [];
                segBelow = []; segAbove = [];
                planesBetween = [];

                s = 0;
                for kk = 1:L
                    % nominal kk
                    s0 = s; s1 = s + h_eff(kk);
                    segS0(end+1,1) = s0; %#ok<AGROW>
                    segS1(end+1,1) = s1;
                    segType(end+1,1) = 0;
                    segK(end+1,1) = kk;
                    segBelow(end+1,1) = 0;
                    segAbove(end+1,1) = 0;

                    s = s1;

                    % interface between kk and kk+1
                    if add_interface && kk < L
                        planesBetween(end+1,1) = s; %#ok<AGROW>  % nominal->interface face

                        s0 = s; s1 = s + i_th;
                        segS0(end+1,1) = s0;
                        segS1(end+1,1) = s1;
                        segType(end+1,1) = 1;
                        segK(end+1,1) = 0;
                        segBelow(end+1,1) = kk;
                        segAbove(end+1,1) = kk + 1;

                        s = s1;

                        planesBetween(end+1,1) = s; %#ok<AGROW>  % interface->nominal face
                    end
                end

                % exclude global ends (not "between layers")
                if ~isempty(planesBetween)
                    planesBetween = planesBetween(planesBetween > tol & planesBetween < (totalH - tol));
                end
            end

            % build stacks
            % Pillar: use effective layers (layer 1 trimmed by int_th/2 for
            % the top interface with ground) to match the actual mesh geometry.
            pillar_layers_eff = obj.pillar_layers(:);
            if isscalar(i_th) && i_th > 0 && ~isempty(obj.ground_layers_eff)
                pillar_layers_eff(1) = pillar_layers_eff(1) - i_th/2;
            end
            [pS0,pS1,pType,pK,pBelow,pAbove,pPlanes,Hp] = buildStack(pillar_layers_eff);

            % Ground: layers go from z=-depth upward, so layer 1 is at BOTTOM (most negative z)
            % and layer L is at TOP (closest to z=0). We use sg = -z so sg=0 at z=0.
            % To make sg=0 correspond to the TOP layer (L), we flip the arrays.
            ground_layers_flipped = flip(obj.ground_layers_eff(:));
            ground_chem_flipped = flip(obj.ground_chem(:));
            [gS0,gS1,gType,gK,gBelow,gAbove,gPlanes,Dg] = buildStack(ground_layers_flipped);

            % ------------------------------------------------------------
            % Special interface slab between ground top layer and pillar layer 1
            % Interface is CENTERED at z=0, spanning z in [-i_th/2, +i_th/2]
            % This matches the rule: interface takes i_th/2 from each adjacent nominal layer
            % ------------------------------------------------------------
            hasTopInterface = isscalar(i_th) && i_th > 0 && ~isempty(obj.pillar_layers) && ~isempty(obj.ground_layers_eff);

            if hasTopInterface
                z_top = i_th/2;      % top of interface (boundary with pillar layer 1)
                z_bot = -i_th/2;     % bottom of interface (boundary with ground top layer)

                inTopIf = (z_local >= (z_bot - tol)) & (z_local <= (z_top + tol));

                if any(inTopIf(:))
                    % Chemistry interpolated from ground top layer (at z=-i_th/2) to pillar layer 1 (at z=+i_th/2)
                    % ground_chem(end) is the topmost ground layer (closest to z=0)
                    cG = obj.ground_chem(end);  % topmost ground layer
                    cP = obj.pillar_chem(1);    % bottom pillar layer

                    if outIsCell
                        % cannot interpolate cells; pick nearest side
                        for ii = find(inTopIf(:))'
                            zi = z_local(ii);
                            u = (zi - z_bot) / max(1e-15, (z_top - z_bot)); % 0 at bottom, 1 at top
                            if u >= 0.5
                                chem{ii} = cP;
                            else
                                chem{ii} = cG;
                            end
                        end
                    else
                        u = (z_local(inTopIf) - z_bot) ./ max(1e-15, (z_top - z_bot)); % 0..1
                        u = max(0, min(1, u));
                        chem(inTopIf) = (1-u).*cG + u.*cP;
                    end

                    % lays_between_layers: true on BOTH boundary faces of top interface
                    % z=+i_th/2 is boundary between top interface and pillar layer 1
                    % z=-i_th/2 is boundary between top interface and ground top layer
                    zTopIf = z_local(inTopIf);
                    onTopFace = abs(zTopIf - z_top) <= tol;
                    onBotFace = abs(zTopIf - z_bot) <= tol;

                    lays_between_layers(inTopIf) = onTopFace | onBotFace;
                end
            end

            % masks for the remaining points (excluding the special top-interface slab)
            doneMask = false(sz);
            if hasTopInterface
                % Top interface is centered at z=0, spanning [-i_th/2, +i_th/2]
                doneMask = doneMask | ((z_local >= (-i_th/2 - tol)) & (z_local <= (i_th/2 + tol)));
            end

            % Pillar: z > i_th/2 (strictly above top interface)
            % Ground: everything else not in top interface
            pilMask = (z_local > (i_th/2 + tol)) & ~doneMask;
            grdMask = ~pilMask & ~doneMask;

            % =========================
            % PILLAR (z > i_th/2)
            % =========================
            if any(pilMask(:))
                zp = z_local(pilMask);
                % Pillar mesh starts at z = +i_th/2, stack coordinate s = z - i_th/2
                % s=0 corresponds to z = +i_th/2 (bottom of pillar layer 1)
                sp = zp - i_th/2;
                sp = max(0, min(Hp, sp));  % clamp to valid stack range [0, Hp]

                % Detect boundary planes (strict tolerance)
                if ~isempty(pPlanes)
                    onPlane = false(size(sp));
                    for j = 1:numel(pPlanes)
                        onPlane = onPlane | (abs(sp - pPlanes(j)) <= tol);
                    end
                    lays_between_layers(pilMask) = onPlane;
                end

                if outIsCell
                    tmp = cell(size(sp));
                else
                    tmp = zeros(size(sp));
                end

                for ii = 1:numel(sp)
                    s = sp(ii);
                    j = find((s >= pS0 - tol) & (s <= pS1 + tol), 1, 'first');

                    if isempty(j)
                        k = numel(obj.pillar_layers);
                        if outIsCell, tmp{ii} = obj.pillar_chem{k}; else, tmp(ii) = obj.pillar_chem(k); end
                        continue;
                    end

                    if pType(j) == 0
                        % nominal layer
                        k = pK(j);
                        if outIsCell, tmp{ii} = obj.pillar_chem{k}; else, tmp(ii) = obj.pillar_chem(k); end
                    else
                        % interface slab - interpolate
                        k0 = pBelow(j); k1 = pAbove(j);
                        u  = (s - pS0(j)) / max(1e-15, (pS1(j) - pS0(j)));
                        u  = max(0, min(1, u));
                        if outIsCell
                            if u >= 0.5
                                tmp{ii} = obj.pillar_chem{k1};
                            else
                                tmp{ii} = obj.pillar_chem{k0};
                            end
                        else
                            c0 = obj.pillar_chem(k0);
                            c1 = obj.pillar_chem(k1);
                            tmp(ii) = (1-u)*c0 + u*c1;
                        end
                    end
                end
                chem(pilMask) = tmp;
            end

            % =========================
            % GROUND (z < -i_th/2)
            % =========================
            if any(grdMask(:))
                zg = z_local(grdMask);
                % Ground mesh extends from z = -Dg - i_th/2 to z = -i_th/2
                % Stack coordinate: s=0 at z = -i_th/2 (top of ground), s=Dg at z = -Dg - i_th/2
                sg = -zg - i_th/2;  % s=0 at z=-i_th/2, s=Dg at z=-Dg-i_th/2
                sg = max(0, min(Dg, sg));  % clamp to valid stack range [0, Dg]

                % Detect boundary planes (strict tolerance)
                if ~isempty(gPlanes)
                    onPlane = false(size(sg));
                    for j = 1:numel(gPlanes)
                        onPlane = onPlane | (abs(sg - gPlanes(j)) <= tol);
                    end
                    lays_between_layers(grdMask) = onPlane;
                end

                if outIsCell
                    tmp = cell(size(sg));
                else
                    tmp = zeros(size(sg));
                end

                for ii = 1:numel(sg)
                    s = sg(ii);
                    j = find((s >= gS0 - tol) & (s <= gS1 + tol), 1, 'first');

                    if isempty(j)
                        % fallback to deepest layer (layer 1 in original = last in flipped)
                        k = numel(ground_chem_flipped);
                        if outIsCell, tmp{ii} = ground_chem_flipped{k}; else, tmp(ii) = ground_chem_flipped(k); end
                        continue;
                    end

                    if gType(j) == 0
                        % nominal layer - use flipped chem array
                        k = gK(j);
                        if outIsCell, tmp{ii} = ground_chem_flipped{k}; else, tmp(ii) = ground_chem_flipped(k); end
                    else
                        % interface slab - interpolate using flipped chem
                        kBelow = gBelow(j);  % layer index below (smaller s)
                        kAbove = gAbove(j);  % layer index above (larger s)
                        u  = (s - gS0(j)) / max(1e-15, (gS1(j) - gS0(j)));
                        u  = max(0, min(1, u));
                        if outIsCell
                            if u >= 0.5
                                tmp{ii} = ground_chem_flipped{kAbove};
                            else
                                tmp{ii} = ground_chem_flipped{kBelow};
                            end
                        else
                            cBelow = ground_chem_flipped(kBelow);
                            cAbove = ground_chem_flipped(kAbove);
                            tmp(ii) = (1-u)*cBelow + u*cAbove;
                        end
                    end
                end
                chem(grdMask) = tmp;
            end
        end



        %==================================================================
        % Z corner coords
        %==================================================================
        function computeZNodalCoords(obj)
            z_tolerance=1.0E-9;
            zCornersPoints = [ ...
                obj.mesh.nodes(obj.mesh.elems(:,1),  :); ...
                obj.mesh.nodes(obj.mesh.elems(:,3),  :); ...
                obj.mesh.nodes(obj.mesh.elems(:,7),  :); ...
                obj.mesh.nodes(obj.mesh.elems(:,9),  :); ...
                obj.mesh.nodes(obj.mesh.elems(:,19), :); ...
                obj.mesh.nodes(obj.mesh.elems(:,21), :); ...
                obj.mesh.nodes(obj.mesh.elems(:,25), :); ...
                obj.mesh.nodes(obj.mesh.elems(:,27), :)  ];

            obj.z_corners_coords = sort(unique(round(zCornersPoints(:,3) / z_tolerance) * z_tolerance), 'descend');
            obj.z_coords = sort(unique(round(obj.mesh.nodes(:,3) / z_tolerance) * z_tolerance), 'descend');
        end

        %==================================================================
        % Eksport siatki i danych do pliku FEAP
        % (bez zmian względem Twojej wersji)
        %==================================================================
        function FEAP_Export( obj, filename ) 
            myfile = fopen(filename, "w");
            if myfile < 0
                error('Nie można otworzyć pliku %s do zapisu.', filename);
            end

            feapNum = [1 3 9 7 19 21 27 25 2 6 8 4 20 24 26 22 10 12 18 16 5 23 13 15 11 17 14];

            % tol = 1E-9;
            %  nodes = round(obj.mesh.nodes / tol) * tol;
            nodes = obj.mesh.nodes;

            fprintf(myfile, "feap * * pillar \n  %d %d 0 3 5 27 \n\n", ...
                size(obj.mesh.nodes,1), size(obj.mesh.elems,1));
            fprintf(myfile, "COORdinates\n");
            for k = 1:size(obj.mesh.nodes, 1)
                fprintf(myfile, "%d   0   %.6f   %.6f   %.5f\n", ...
                    k, nodes(k,1), nodes(k,2), nodes(k,3));
            end

            fprintf(myfile, "\n\nELEMents\n");
            for k = 1:size(obj.mesh.elems,1)
                felems = obj.mesh.elems(k, feapNum);
                fprintf(myfile, "%d 0 1", k);
                for e = 1:numel(feapNum)
                    fprintf(myfile, " %d", felems(e));
                    if e == 13
                        fprintf(myfile, "\n");
                    end
                end
                fprintf(myfile, "\n");
            end
            fprintf(myfile, "\n");

            fprintf(myfile, "\n BOUNdary");
            fprintf(myfile, "\n1 1 0 0 0  -1 -1");
            fprintf(myfile, "\n %d 1 0 0 0  1  1\n", size(obj.mesh.nodes,1));

            fprintf(myfile, "\ngap 0.001");
            fprintf(myfile, "\nEBOUndary ADD\n");
            fprintf(myfile, "1 0.0   1 0 0  1 1 ! plane x = 0  -> fix u_x\n");
            fprintf(myfile, "2 0.0   0 1 0  1 1 ! plane y = 0  -> fix u_y\n");
            fprintf(myfile, "3 %7.5E   1 1 1  1 1  ! plane z = min  -> fix u_x,u_y,u_z\n", ...
                min(obj.mesh.nodes(:,3)));

            [chem_from_z, lays_between_layers] = obj.chemFromZ(obj.z_coords);
            fprintf(myfile, "\n EDIS\n");
            for k = 1:size(obj.z_coords,1)
                str_layer="";
                if lays_between_layers(k)
                    str_layer=" ! ------layer surface --------------------------------";
                end
                if obj.z_coords(k)>0
                    fprintf(myfile, "  3   %.5f  0  0  0  %1.2f 0.0   %s\n", ...
                     obj.z_coords(k), chem_from_z(k),str_layer);
                else
                    fprintf(myfile, "  3   %.5f  0  0  0  0.0 %1.2f   %s\n", ...
                     obj.z_coords(k), chem_from_z(k),str_layer);
                end
            end
            fprintf(myfile, "\n");

            % dalsza konfiguracja FEAP (jak u Ciebie)
            fprintf(myfile, "mate,1\n");
            fprintf(myfile, "user,14\n");
            fprintf(myfile, "3,2,2,2,30,0.,1.\n");
            fprintf(myfile, "2 -1 -1 0  0 0 0 1\n");
            fprintf(myfile, "3.189E-10 5.185E-10\n");
            fprintf(myfile, "div(sig) -2 1   1\n");
            fprintf(myfile, "0.\n");
            fprintf(myfile, "390.0d9,145.0d9,106.0d9,398.0d9,105.0d9\n");
            fprintf(myfile, "223.0d9,115.0d9,92.0d9,224.0d9,50.0d9\n");
            fprintf(myfile, "195.0d9, 72.5d9, 53.0d9,199.0d9, 52.5d9\n");
            fprintf(myfile, "div(x_n) -1 1  -4\n");
            fprintf(myfile, "3.533E-10   5.693E-10\n");
            fprintf(myfile, "3.189E-10 5.185E-10\n");
            fprintf(myfile, "\n");
            fprintf(myfile, "end\n\n");
            fprintf(myfile, "TIE\n\n");
            fprintf(myfile, "BATCh\n");
            fprintf(myfile, "PROP\n");
            fprintf(myfile, "DT,,1\n");
            fprintf(myfile, "END\n");
            fprintf(myfile, "2 2\n");
            fprintf(myfile, "0 0 1 1\n\n");
            fprintf(myfile, "batch\n");
            fprintf(myfile, "plot,pers,1\n");
            fprintf(myfile, "plot,hide\n");
            fprintf(myfile, "plot,fill\n");
            fprintf(myfile, "plot,defo\n");
            fprintf(myfile, "plot,mesh\n");
            fprintf(myfile, "plot,load\n");
            fprintf(myfile, "plot,axis\n");
            fprintf(myfile, "plot,cont,4\n");
            fprintf(myfile, "end\n");
            fprintf(myfile, "0\n");
            fprintf(myfile, "-2000. -4000. 2000.\n");
            fprintf(myfile, "  0.  0. 2.\n\n");
            fprintf(myfile, "batch\n");
            fprintf(myfile, "plot,mesh\n");
            fprintf(myfile, "plot,defo,1,1\n");
            fprintf(myfile, "plot,cont,4\n");
            fprintf(myfile, "end\n\n");
            fprintf(myfile, "batch\n");
            fprintf(myfile, "opti\n\n");
            fprintf(myfile, "LOOP,,1\n");
            fprintf(myfile, "TIME\n");
            fprintf(myfile, "LOOP,,99\n");
            fprintf(myfile, "utan,,1\n");
            fprintf(myfile, "plot,cont,4\n");
            fprintf(myfile, "plot,stre,5\n");
            fprintf(myfile, "plot,stre,4\n");
            fprintf(myfile, "plot,stre,3\n");
            fprintf(myfile, "plot,stre,2\n");
            fprintf(myfile, "plot,stre,1\n");
            fprintf(myfile, "!plot,stre,6\n");
            fprintf(myfile, "NEXT\n");
            fprintf(myfile, "save\n");
            fprintf(myfile, "disp,all\n");
            fprintf(myfile, "stre,all\n");
            fprintf(myfile, "stre,node,all\n");
            fprintf(myfile, "NEXT\n\n");
            fprintf(myfile, "end\n\n");
            fprintf(myfile, "inte\n");
            fprintf(myfile, "stop\n");

            fclose(myfile);
        end

        function report = checkMeshIntegrity(obj, tol)
            % checkMeshIntegrity  Basic integrity checks for hex meshes (L27 supported).
            %
            % report fields:
            %   .nNodes, .nElems
            %   .hasNaNInf
            %   .badElemIndexCount
            %   .unusedNodeCount
            %   .duplicateNodePairsCount
            %   .zeroCornerEdgeCount
            %   .degenerateElemCount
            %   .invertedElemCount
            %   .nComponents
            %   .largestComponentFrac
            %
            % Usage:
            %   r = model.checkMeshIntegrity(1e-6);
            %   disp(r)
    
            if nargin < 2 || isempty(tol)
                tol = 1e-9;
            end
    
            X = obj.mesh.nodes;
            E = obj.mesh.elems;
    
            report = struct();
            report.nNodes = size(X,1);
            report.nElems = size(E,1);
    
            % ---------------------------------
            % 0) NaN/Inf coordinates
            % ---------------------------------
            report.hasNaNInf = any(~isfinite(X(:)));
    
            % ---------------------------------
            % 1) Element connectivity index validity
            % ---------------------------------
            minIdx = min(E(:));
            maxIdx = max(E(:));
            badIdxMask = (E(:) < 1) | (E(:) > report.nNodes) | ~isfinite(E(:));
            report.badElemIndexCount = nnz(badIdxMask);
            report.elemIndexRange = [minIdx, maxIdx];
    
            % ---------------------------------
            % 2) Unused nodes
            % ---------------------------------
            used = false(report.nNodes,1);
            if report.badElemIndexCount == 0
                used(unique(E(:))) = true;
            end
            report.unusedNodeCount = nnz(~used);
    
            % ---------------------------------
            % 3) Duplicate nodes (within tolerance)
            %    Uses rounding grid; counts duplicates (approx).
            % ---------------------------------
            t = tol;
            if isprop(obj,'z_tolerance') && ~isempty(obj.z_tolerance)
                % Honor the caller tolerance as the primary diagnostic setting.
                % z_tolerance is a lower-bound for snapping/classification,
                % not a request to silently tighten duplicate reporting.
                t = max(t, obj.z_tolerance);
            end
            key = round(X ./ t);
            [~, ia, ic] = unique(key, 'rows', 'stable');
            report.duplicateNodePairsCount = report.nNodes - numel(ia);
    
            % ---------------------------------
            % 4) Degenerate / inverted elements (corner-based signed volume)
            %    Works for L27 as long as corner indices are [1 3 7 9 19 21 25 27]
            % ---------------------------------
            if size(E,2) < 8
                report.degenerateElemCount = NaN;
                report.invertedElemCount   = NaN;
            else
                if size(E,2) >= 27
                    c = [1 3 7 9 19 21 25 27];     % your L27 corner pattern
                else
                    c = 1:8;                       % fallback
                end
    
                Ec = E(:,c);
                % tetra split of hex corners (gives signed volume)
                tets = [1 2 4 5;
                        2 3 4 7;
                        2 4 7 5;
                        2 7 6 5;
                        4 7 8 5];
    
                vols = nan(report.nElems,1);
    
                for e = 1:report.nElems
                    idx = Ec(e,:);
                    if any(idx < 1) || any(idx > report.nNodes)
                        continue;
                    end
                    xe = X(idx,:);
    
                    v = 0;
                    for k = 1:size(tets,1)
                        a = xe(tets(k,1),:);
                        b = xe(tets(k,2),:);
                        c2= xe(tets(k,3),:);
                        d = xe(tets(k,4),:);
                        v = v + det([b-a; c2-a; d-a]) / 6;
                    end
                    vols(e) = v;
                end

                scaleX = max(1, max(abs(X(:))));
                volTol = 1e-14 * scaleX;

                finiteVols = isfinite(vols);
                strong = finiteVols & (abs(vols) > volTol);

                % Normalize global sign convention: if most strong volumes are
                % negative, flip all signs before counting inversions.
                if nnz(strong) > 0
                    if nnz(vols(strong) < 0) > nnz(vols(strong) > 0)
                        vols = -vols;
                    end
                end

                bad = ~isfinite(vols) | (abs(vols) <= volTol);
                inv = vols < -volTol;

                report.degenerateElemCount = nnz(bad);
                report.invertedElemCount   = nnz(inv);
                report.minSignedVolume     = min(vols);
                report.maxSignedVolume     = max(vols);
            end
    
            % ---------------------------------
            % 5) Zero-length corner edges (quick indicator of collapse)
            % ---------------------------------
            report.zeroCornerEdgeCount = 0;
            if size(E,2) >= 8 && report.badElemIndexCount == 0
                if size(E,2) >= 27
                    c = [1 3 7 9 19 21 25 27];
                else
                    c = 1:8;
                end
                Ec = E(:,c);
    
                % corner edges of a hex (in corner ordering)
                edges = [1 2; 2 3; 3 4; 4 1;   % bottom loop
                         5 6; 6 7; 7 8; 8 5;   % top loop
                         1 5; 2 6; 3 7; 4 8];  % verticals
    
                Xc = X; %#ok<NASGU>
                zc = 0;
                for e = 1:report.nElems
                    idx = Ec(e,:);
                    xe = X(idx,:);
                    for k = 1:size(edges,1)
                        a = xe(edges(k,1),:);
                        b = xe(edges(k,2),:);
                        if norm(a-b) <= tol
                            zc = zc + 1;
                        end
                    end
                end
                report.zeroCornerEdgeCount = zc;
            end
    
            % ---------------------------------
            % 6) Connected components (node adjacency via elements)
            % ---------------------------------
            if report.badElemIndexCount == 0
                n = report.nNodes;
                adj = cell(n,1);
                for e = 1:report.nElems
                    en = unique(E(e,:));
                    en = en(en >= 1 & en <= n);
                    for ii = 1:numel(en)
                        ni = en(ii);
                        adj{ni} = [adj{ni}, en]; %#ok<AGROW>
                    end
                end
                for i = 1:n
                    adj{i} = unique(adj{i}(adj{i} ~= i));
                end
    
                comp = zeros(n,1);
                cid = 0;
                for i = 1:n
                    if ~used(i) || comp(i) ~= 0
                        continue
                    end
                    cid = cid + 1;
                    q = i;
                    comp(i) = cid;
                    while ~isempty(q)
                        u = q(end); q(end) = [];
                        nb = adj{u};
                        nb = nb(used(nb));
                        nb = nb(comp(nb) == 0);
                        comp(nb) = cid;
                        q = [q; nb(:)]; %#ok<AGROW>
                    end
                end
    
                report.nComponents = cid;
                if cid == 0
                    report.largestComponentFrac = 0;
                else
                    counts = accumarray(comp(comp>0), 1);
                    report.largestComponentFrac = max(counts) / nnz(used);
                end
            else
                report.nComponents = NaN;
                report.largestComponentFrac = NaN;
            end
    
            % ---------------------------------
            % Summary flag
            % ---------------------------------
            report.ok = ...
                ~report.hasNaNInf && ...
                report.badElemIndexCount == 0 && ...
                report.invertedElemCount == 0 && ...
                report.degenerateElemCount == 0 && ...
                (isnan(report.nComponents) || report.nComponents <= 1);

            [dmin, pair] = obj.mesh.minNodeDistance();

    
            % Print concise summary (optional)
            fprintf("Mesh integrity:\n");
            fprintf("  Nodes: %d, Elems: %d\n", report.nNodes, report.nElems);
            fprintf("  NaN/Inf coords: %d\n", report.hasNaNInf);
            fprintf("  Bad elem indices: %d (range [%g..%g])\n", report.badElemIndexCount, report.elemIndexRange(1), report.elemIndexRange(2));
            fprintf("  Unused nodes: %d\n", report.unusedNodeCount);
            fprintf("  Duplicate nodes (approx): %d\n", report.duplicateNodePairsCount);
            fprintf("  Degenerate elems: %d, Inverted elems: %d\n", report.degenerateElemCount, report.invertedElemCount);
            fprintf("  Components: %g, Largest component frac: %g\n", report.nComponents, report.largestComponentFrac);
            fprintf("  Minimal distance: %g between nodes: %d, %d \n", dmin, pair(1), pair(2));      
            obj.mesh.nodes(pair(2),:)
            obj.mesh.nodes(pair(1),:)
            fprintf("  Zero-length corner edges: %d\n", report.zeroCornerEdgeCount);
            fprintf("  OK: %d\n", report.ok);
        end

        function plotBadHexFaces(obj, badE, faceAlpha)
            if nargin < 3, faceAlpha = 0.25; end
            E = obj.mesh.elems(badE,:);
            X = obj.mesh.nodes;
        
            c = [1 3 7 9 19 21 25 27];   % L27 corners
            Ec = E(:,c);
        
            % faces in terms of the 8 corners (hex)
            F = [1 2 3 4;
                 5 6 7 8;
                 1 2 6 5;
                 2 3 7 6;
                 3 4 8 7;
                 4 1 5 8];
        
            faces = zeros(size(Ec,1)*6, 4);
            for i = 1:6
                faces( (i-1)*size(Ec,1)+ (1:size(Ec,1)), : ) = Ec(:,F(i,:));
            end
        
            patch('Vertices', X, 'Faces', faces, ...
                  'FaceAlpha', faceAlpha, 'EdgeAlpha', 0.4);
        end

        function nodes = snapCriticalZCoordinates(obj, nodes)
            % Snap nodes to exact z-coordinates at layer boundaries
            tol = obj.z_tolerance;
            
            % Critical z-values
            ground_depth = sum(obj.ground_layers_eff);
            pillar_height = sum(obj.pillar_layers);
            
            critical_z = [
                -ground_depth;  % bottom of ground
                -obj.int_th/2;                  % bottom of interface
                0.0;                            % center (NOT a boundary)
                obj.int_th/2;                   % top of interface
                pillar_height;  % top of pillar
            ];
            
            % Add layer boundaries
            z_ground = -obj.int_th/2;
            for i = numel(obj.ground_layers_eff):-1:1
                z_ground = z_ground -  obj.ground_layers_eff(i);
                critical_z = [critical_z; z_ground]; %#ok<AGROW>
            end
            
            z_pillar = obj.int_th/2;
            for i = 1:numel(obj.pillar_layers)
                z_pillar = z_pillar + obj.pillar_layers(i);
        critical_z = [critical_z; z_pillar]; %#ok<AGROW>
    end
    
    critical_z = unique(critical_z);
    
    % Snap nodes near critical z-values
    z = nodes(:,3);
    for i = 1:numel(critical_z)
        mask = abs(z - critical_z(i)) < tol;
        nodes(mask, 3) = critical_z(i);
    end
end

        function nodes = snapCriticalRCoordinates(obj, nodes, critical_radii, tol)
            if nargin < 4, tol = obj.z_tolerance; end
            
            x = nodes(:,1);
            y = nodes(:,2);
            r = hypot(x, y);
            
            for i = 1:numel(critical_radii)
                R_crit = critical_radii(i);
                mask = abs(r - R_crit) < tol;
                if any(mask)
                    % Scale x,y to exact radius
                    r_mask = r(mask);
                    scale = R_crit ./ r_mask;
                    nodes(mask, 1) = x(mask) .* scale;
                    nodes(mask, 2) = y(mask) .* scale;
                end
            end
        end

        function nodes = snapSymmetryPlanes(obj, nodes)
            tol = obj.z_tolerance;
            
            % Snap to x=0 plane (theta = 90°)
            mask_x = abs(nodes(:,1)) < tol;
            nodes(mask_x, 1) = 0.0;
            
            % Snap to y=0 plane (theta = 0°)
            mask_y = abs(nodes(:,2)) < tol;
            nodes(mask_y, 2) = 0.0;
        end

         function obj = snapAllCoordinates(obj)
            % Snap all coordinates to ensure conformity
            
            % 1. Z-coordinates
            obj.mesh.nodes = obj.snapCriticalZCoordinates(obj.mesh.nodes);
            
            % 2. Symmetry planes
            obj.mesh.nodes = obj.snapSymmetryPlanes(obj.mesh.nodes);
            
            % 3. Critical radii
            pillar_height = sum(obj.pillar_layers);
            Rbase = obj.top_R + pillar_height * tan(5*pi/180);
            Rout = 1.2*Rbase;
            Rbank = 1.5*Rout;
            Rtile = 1.2*Rbank;
            
            critical_radii = [obj.top_R, Rbase, Rout, Rbank, Rtile];
            obj.mesh.nodes = obj.snapCriticalRCoordinates(obj.mesh.nodes, critical_radii, obj.z_tolerance);
        end
    end

    methods(Static)

       

        function X = snapR(X, R, tol)
            x = X(:,1); y = X(:,2);
            r = hypot(x,y);
            m = abs(r - R) <= tol & r > 0;
            s = R ./ r(m);
            x(m) = x(m).*s;
            y(m) = y(m).*s;
            X(:,1)=x; X(:,2)=y;
        end

        function reportSeamNearRadius(mesh, r0, tol)
            X = mesh.nodes;
            r = hypot(X(:,1), X(:,2));
            idx = find(abs(r - r0) < tol);
            XY = X(idx,1:2);
            key = round(XY / tol);
            [~,~,ic] = unique(key,'rows');
            counts = accumarray(ic,1);
            fprintf("r=%.6g: nodes=%d, clustered=%d, duplicates=%d\n", ...
                r0, numel(idx), numel(counts), nnz(counts>1));
        end
        
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
            % Apply smooth "bowl" depression on the TOP surface of a ring.
            %
            % IMPORTANT:
            %   We only move nodes on z == z_top. This preserves flat,
            %   horizontal internal interfaces/layers below, which is
            %   required by the layered-process constraint.
            
                x = nodes(:,1) - x0_xy(1);
                y = nodes(:,2) - x0_xy(2);
                z = nodes(:,3);
            
                r = hypot(x,y);
            
                % --- tolerances (scale-aware) ---
                epsR = max(1e-10*max(1,Rout), 1e-8);              % radial tol
                epsZ = max(1e-10*max(1,abs(z_top-z_bot)), 1e-8);  % z tol
            
                inR = (r >= Rin - epsR) & (r <= Rout + epsR);
                onTop = abs(z - z_top) <= epsZ;
            
                if ~any(inR & onTop)
                    return;
                end
            
                onInner = abs(r - Rin)  <= epsR;
                onOuter = abs(r - Rout) <= epsR;
            
                % Only top-surface interior of the ring
                mask = inR & onTop & ~(onInner | onOuter);
            
                if ~any(mask)
                    return;
                end
            
                t = (r(mask) - Rin) ./ max(1e-15, (Rout - Rin));
                t = max(0, min(1, t));
            
                % radial bowl shape with zero slope at boundaries:
                % 0 at t=0 and t=1, max at t=0.5, C1-continuous at edges.
                bowl = sin(pi .* t).^2;
            
                dz = -depth .* bowl;
            
                % apply
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
		            Rin,  thB, z1;
                    RoutB,thB, z1;

                    Rin,  thA, z2;
                    RoutA,thA, z2;
                    Rin,  thB, z2;
                    RoutB,thB, z2;
                ];
        
                m = Mesh();
                % NOTE: adjust argument order to match YOUR Mesh.addRing3D signature.
                % Your earlier call suggests: addRing3D(x0, obj_sf, x, nr, nz, nc, localNodes)
                m.addRing3D([0 0 0], obj_sf, x, nr, 1, nz, sf.localNodes);
        
                mesh.mergeMesh(m);
            end
        end


    end
end

classdef ManipulatorModel3D < handle
 
    
    properties
        analysis, frame_analysis, fixedEdgeSelector, alpha, R, r, mesh, frame_mesh, elems, fe, xEnd, frameNodes, frameElem
        const_elems, upper_nodes, loadSurfaceNodes, fixedSurfaceNodes, halfSegmentNelems, use_offset, qnodal_solid, qnodal_top;
        constEndRing, constMiddleRing;
        couplingSections, debugCoupling;
        % Mesh resolution stored for linkage diagnostics and rotation-aware mode.
        resCirc, resLen, resTh, nCircDiv;
        % Fillet radius applied at each interior bevel joint (0 = no fillet).
        filletR;
    end

    methods
        function obj = ManipulatorModel3D(E,nu,ls,R,r, res, res_th, alpha, betas, ShapeFn, use_offset, Pz, constEndRing, constMiddleRing, nCircDiv, filletR)
            if nargin < 13 || isempty(constEndRing)
                constEndRing = true;
            end
            if nargin < 14 || isempty(constMiddleRing)
                constMiddleRing = true;
            end
            if nargin < 15 || isempty(nCircDiv)
                nCircDiv = 1;
            end
            if nargin < 16 || isempty(filletR)
                filletR = 0;
            end
            obj.nCircDiv = nCircDiv;
            obj.filletR  = filletR;
            alpha=alpha*pi/180;
            obj.alpha = alpha;
            obj.R = R;
            obj.r = r;
            obj.constEndRing = logical(constEndRing);
            obj.constMiddleRing = logical(constMiddleRing);
            obj.debugCoupling = false;
            betas=betas*pi/180;

            if use_offset
                obj.use_offset = 1;
            else
                obj.use_offset = 0;
            end

            obj.mesh=Mesh();
            obj.geterateManipulator( ls, R, r, res, res_th, alpha, betas, ShapeFn, nCircDiv, filletR );
            obj.fe = SolidElasticElem( ShapeFn, obj.elems );

            obj.fe.props.h=1;
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material);

            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
           
            fixedEdgeSelector = Selector( obj.fixedSurfaceNodes );
            loadedFaceSelector = Selector( obj.loadSurfaceNodes );
            A = pi*(R^2 - r^2);
            
            % Keep the direct solid reference load aligned with the frame
            % tip load so coupled-vs-solid comparisons use the same forcing.
            obj.analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 -Pz/A] ));
            % Clamp the full root section in the solid model so the
            % reference 3D problem is mechanically equivalent to the frame
            % root clamp. This also removes the artificial single-node anchor.
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"], [0 0 0] );

            frame_elems = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
            obj.frame_mesh=Mesh();
            obj.frame_mesh.nodes=obj.frameNodes;
            obj.frame_mesh.elems = frame_elems;
            
            obj.frameElem=Frame3D(frame_elems,E,nu,R,r);
            obj.frame_analysis = LinearElasticityWeighted( obj.frameElem, obj.frame_mesh, false );
            obj.frame_analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 0 0 0 0]);
            obj.frame_analysis.loadClosestNode(obj.frame_mesh.nodes(end,:), ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 -Pz 0 0 0] );

            % Optional debug hook: print the resolved frame DOF order so
            % coupling diagnostics can confirm the frame layout explicitly.
            if obj.debugCoupling
                disp('frame_analysis.ndofs =');
                disp(obj.frame_analysis.ndofs);
            end

            % Temporary debug hook: verify the frame tip load acts only on uz.
            tipFrameNode = obj.frame_mesh.findClosestNode(obj.frame_mesh.nodes(end,:));
            iuzFrame = obj.frame_analysis.findDOFsIndices("uz");
            tipLoad = obj.frame_analysis.Pnodal(tipFrameNode,:);
            otherLoadMask = true(1, numel(obj.frame_analysis.ndofs));
            otherLoadMask(iuzFrame) = false;
            assert(abs(tipLoad(iuzFrame) + Pz) < 1.0e-12, ...
                'Frame tip load must be applied on uz only.');
            assert(norm(tipLoad(otherLoadMask)) < 1.0e-12, ...
                'Unexpected non-uz frame tip load detected.');

        end

        function ns = determineSegment( obj, elem )
            nhs = round(elem ./ obj.halfSegmentNelems );  
            ns = ceil((nhs - 1) / 2) + 1;
        end

        function geterateManipulator(obj, ls, R, r, res, res_th, alpha, betas, sf, nCircDiv, filletR)
            if nargin < 10 || isempty(nCircDiv)
                nCircDiv = 1;
            end
            if nargin < 11 || isempty(filletR)
                filletR = 0;
            end
            Th = R - r;
            obj.couplingSections = struct('frameNode', {}, 'adjacentFrameNode', {}, ...
                'nodes', {}, 'label', {});

            % --- thickness resolution (radial) ---
            resTh = max(1, round(res_th));

            % --- choose circum + length resolution consistently with resTh ---
            % (your original logic: scale with Th so elements stay ~ isotropic)
            resCirc = max(3, round(2*pi*R/Th * resTh));
            resLen  = max(1, round(ls/Th * resTh));

            % --- rotation-aware mode (use_offset=false): snap resCirc so every
            %     joint angle is an exact integer multiple of 2*pi/resCirc.
            %     Junction nodes then align without the phase trick, and element e
            %     covers the same local circumferential band in every configuration.
            %     nCircDiv must be > 1 (use arm.nCircDiv=8 for the standard configs).
            if ~obj.use_offset
                assert(nCircDiv > 1, ...
                    ['ManipulatorModel3D: use_offset=false requires nCircDiv > 1 so that ' ...
                     'resCirc is snapped and junction nodes align. Pass arm.nCircDiv (e.g. 8).']);
                resCirc = round(resCirc / nCircDiv) * nCircDiv;
                resCirc = max(resCirc, nCircDiv);
            end

            obj.resCirc = resCirc;
            obj.resLen  = resLen;
            obj.resTh   = resTh;
        
            % rotations...
            c=cos(alpha); s=sin(alpha);
            rotCut=[c 0 s; 0 1 0; -s 0 c]';
            c=cos(2*alpha); s=sin(2*alpha);
            rotCut2=[c 0 s; 0 1 0; -s 0 c]';
            c=cos(betas(1)); s=sin(betas(1));
            rotBeta=[c -s 0; s c 0; 0 0 1]';
        
            phase = -betas(1)*obj.use_offset;
        
            % --- first segment ---
            [obj.mesh, obj.elems, startSection1, endSection1] = ...
                obj.generateSegment2a(R, r, ls, phase, rotCut, sf, resTh, resCirc, resLen);
        
            obj.halfSegmentNelems = size(obj.elems,1);
        
            % elements in ONE end slice (k = last along length):
            sliceElems = resTh * resCirc;
            ne = size(obj.elems,1);
            obj.const_elems = obj.localConstRingIds(1, ne, sliceElems);
        
            obj.mesh.nodes = obj.mesh.nodes * rotBeta;
            selector = Selector(@(x)(x(:,3) < 1.0E-4));
            obj.fixedSurfaceNodes = selector.select(obj.mesh.nodes);
        
            prevRot = rotCut * rotBeta;
            obj.xEnd = [0 0 ls];
            xEnds = [[0 0 0]; obj.xEnd];

            % The first half-segment contributes the root section and the
            % first coupling section at frame node 2.
            obj.registerCouplingSection(1, 2, startSection1, "root_section");
            obj.registerCouplingSection(2, 1, endSection1, "joint_2_in");
        
            for k = 2:length(betas)
                c=cos(betas(k)); s=sin(betas(k));
                rotBeta=[c -s 0; s c 0; 0 0 1]';
        
                phase = obj.use_offset * (phase - betas(k));
        
                [mesh1, elems1, startSectionLocal, endSectionLocal] = ...
                    obj.generateSegment2b(R, r, ls, phase, rotCut, sf, resTh, resCirc, resLen);
                [mesh2, elems2, ~, endSection2] = ...
                    obj.generateSegment2a(R, r, ls, phase, rotCut, sf, resTh, resCirc, resLen);
        
                if k < length(betas)
                    mergedElems2 = mesh1.merge(mesh2.nodes, elems2);
                    endSectionLocal = obj.mapNodesByCoordinates(mesh1.nodes, mesh2.nodes(endSection2, :), mesh1.tolerance);
                    elems1 = [elems2; mergedElems2];
                end
        
                mesh1.nodes = (mesh1.nodes + [0 0 ls]) * rotCut * rotBeta * prevRot + obj.xEnd;
                startSectionCoords = mesh1.nodes(startSectionLocal, :);
                endSectionCoords = mesh1.nodes(endSectionLocal, :);
                prevElemCount = size(obj.elems, 1);
                obj.elems = [obj.elems; obj.mesh.merge(mesh1.nodes, elems1)];

                nextFrameNode = size(xEnds, 1) + 1;
                startSectionGlobal = obj.mapNodesByCoordinates(obj.mesh.nodes, startSectionCoords, obj.mesh.tolerance);
                endSectionGlobal = obj.mapNodesByCoordinates(obj.mesh.nodes, endSectionCoords, obj.mesh.tolerance);
                obj.registerCouplingSection(k, nextFrameNode, startSectionGlobal, "joint_" + string(k) + "_out");
                obj.registerCouplingSection(nextFrameNode, k, endSectionGlobal, "joint_" + string(nextFrameNode) + "_in");
        
                % update const_elems for this newly appended part:
                ne = size(obj.elems,1);
                for blockStart = prevElemCount + 1:obj.halfSegmentNelems:ne
                    blockEnd = min(blockStart + obj.halfSegmentNelems - 1, ne);
                    obj.const_elems = [obj.const_elems; ...
                        obj.localConstRingIds(blockStart, blockEnd, sliceElems)]; %#ok<AGROW>
                end
        
                if k < length(betas)
                    obj.xEnd = obj.xEnd + [0 0 2*ls] * rotCut * rotBeta * prevRot;
                else
                    obj.xEnd = obj.xEnd + [0 0 ls]   * rotCut * rotBeta * prevRot;
                end
                xEnds = [xEnds; obj.xEnd];
                prevRot = rotCut2 * rotBeta * prevRot;
            end
        
            tNodes = (obj.mesh.nodes - repmat(obj.xEnd, size(obj.mesh.nodes,1), 1)) * prevRot' * rotCut;
            obj.loadSurfaceNodes = abs(tNodes(:,3)) < 1.0E-4;
            obj.frameNodes = xEnds;
            nCopies = size(obj.elems, 1) / obj.halfSegmentNelems;
            obj.const_elems = obj.fullArmPhysicalRingIds(round(nCopies), sliceElems);

            obj.applyJointFillets(filletR, xEnds, R);
        end

        function ids = localConstRingIds(obj, firstElem, lastElem, sliceElems)
            ids = zeros(0, 1);
            if obj.constEndRing
                ids = [ids; (firstElem:firstElem + sliceElems - 1)']; %#ok<AGROW>
            end
            if obj.constMiddleRing
                ids = [ids; (lastElem - sliceElems + 1:lastElem)']; %#ok<AGROW>
            end
        end

        function ids = fullArmPhysicalRingIds(obj, nCopies, sliceElems)
            H = obj.halfSegmentNelems;
            ids = zeros(0, 1);

            if obj.constEndRing
                ids = [ids; (1:sliceElems)']; %#ok<AGROW>
                for b = 2:2:nCopies
                    blockStart = (b - 1) * H + 1;
                    blockEnd = blockStart + H - 1;
                    ids = [ids; (blockEnd - sliceElems + 1:blockEnd)']; %#ok<AGROW>
                end
            end

            if obj.constMiddleRing
                for b = 1:2:nCopies-1
                    blockStart = (b - 1) * H + 1;
                    blockEnd = blockStart + H - 1;
                    ids = [ids; (blockEnd - sliceElems + 1:blockEnd)']; %#ok<AGROW>
                end
            end

            ids = unique(ids(:));
        end


        function [mesh, elems, startSectionNodes, endSectionNodes] = generateSegment2a(obj, R, r, ls, phase, Redge, sf, resTh, resCirc, resLen)
            Th = R - r;
        
            mesh = Mesh();
            elems = mesh.addRectMesh3D( ...
                R - Th, phase, 0, ...
                Th, 2*pi, 1, ...
                resTh, resCirc, resLen, sf.localNodes);
        
            elems = mesh.transformToCylindrical3D([0 0]);
            axialParam = mesh.nodes(:,3);
            snapTol = 1.0 / mesh.tolerance;
            startSectionNodes = int32(find(abs(axialParam - min(axialParam)) <= snapTol));
            endSectionNodes = int32(find(abs(axialParam - max(axialParam)) <= snapTol));
        
            nodesR1 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)] * Redge;
            nodes   = [nodesR1(:,1), nodesR1(:,2), mesh.nodes(:,3) .* (nodesR1(:,3) + ls)];
            mesh.nodes = nodes;
            elems = mesh.elems;
        end


        function [mesh, elems, startSectionNodes, endSectionNodes] = generateSegment2b(obj, R, r, ls, phase, Redge, sf, resTh, resCirc, resLen)
            Th = R - r;
        
            mesh = Mesh();
            elems = mesh.addRectMesh3D( ...
                R, phase, 0, ...
                -Th, -2*pi, 1, ...
                resTh, resCirc, resLen, sf.localNodes);
        
            elems = mesh.transformToCylindrical3D([0 0]);
            axialParam = mesh.nodes(:,3);
            snapTol = 1.0 / mesh.tolerance;
            startSectionNodes = int32(find(abs(axialParam - min(axialParam)) <= snapTol));
            endSectionNodes = int32(find(abs(axialParam - max(axialParam)) <= snapTol));
        
            nodesR1 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)] * Redge';
            nodes   = [nodesR1(:,1), nodesR1(:,2), (1 - mesh.nodes(:,3)) .* (nodesR1(:,3) - ls)];
            mesh.nodes = nodes;
            elems = mesh.elems;
        end

        function xnew=rotatePoints(x0,R,x)
            xnew=(x-x0)*R+x0;
        end 

        function compute(obj, x)
            if size(x,1)==1
                x=ones(obj.analysis.getTotalElemsNumber(),1);
            end
            obj.analysis.solveWeighted(x);
            obj.analysis.computeElementResults(x);
        end

        function plot(obj)      
            obj.analysis.felems{1}.plot(obj.mesh.nodes);
            % obj.analysis.plotCurrentLoad();
            obj.analysis.plotSupport();
        end

        function plotMesh(obj)       
            obj.fe.face_alpha = 0.2;
            ec = obj.fe.edge_color;
            obj.fe.edge_color='none';
            obj.fe.face_color=[0.4 0.4 0.4];


            obj.fe.plot(obj.mesh.nodes);
            obj.frameElem.plot(obj.frame_mesh.nodes);
            
            obj.fe.face_alpha = 1.0;
            obj.fe.edge_color=ec;
        end

        function maxHM = plotHM_map(obj, x)
            if size(x,1)==1
                x=ones(obj.analysis.getTotalElemsNumber(),1);
            end
            selems = x > 0.5; 
            maxHM = max(obj.analysis.felems{1}.results.gp.all(13,selems,:));
            figure
            hold on, axis on; 
            daspect([1 1 1]);
            obj.analysis.felems{1}.selectedElems=selems;
            obj.analysis.plotMaps(["sHM"],0.0)
            obj.analysis.felems{1}.selectedElems=[];
        end

         function FN = computeFrameNodesOnly(obj,  ls, alpha_deg, betas_deg)
            alpha = deg2rad(alpha_deg);
            betas = deg2rad(betas_deg(:).');   % row
            nJ = numel(betas);
    
            % rotCutT = Ry(alpha)'; rotCut2T = Ry(2*alpha)';
            ca = cos(alpha);  sa = sin(alpha);
            rotCutT  = [ ca  0  sa;  0  1  0; -sa  0  ca ]';
            c2 = cos(2*alpha); s2 = sin(2*alpha);
            rotCut2T = [ c2  0  s2;  0  1  0; -s2  0  c2 ]';
    
            % Rz1T
            cb = cos(betas(1)); sb = sin(betas(1));
            Rz1T = [ cb -sb 0; sb cb 0; 0 0 1 ]';
    
            prevRot = rotCutT * Rz1T;
    
            xEnd = [0 0 ls];
            FN = zeros(nJ,3);
            FN(1,:) = [0 0 0];
            FN(2,:) = xEnd;
    
            for k = 2:nJ
                cb = cos(betas(k)); sb = sin(betas(k));
                RzTk = [ cb -sb 0; sb cb 0; 0 0 1 ]';
                step = (k < nJ) * 2*ls + (k == nJ) * ls;
                xEnd = xEnd + [0 0 step] * rotCutT * RzTk * prevRot;
                FN(k,:) = xEnd;
                % prevRot = Ry(2*alpha)' * Rz(betas(k))' * prevRot
                prevRot = rotCut2T * RzTk * prevRot;
            end
        end
    
        function FNb = computeFrameNodesBatch(obj, ls, alpha_deg, betas_deg_batch)
            B = betas_deg_batch;
            if size(B,1) < size(B,2) && size(B,1) <= 8
                B = B.';
            end
            [nSamples, nJoints] = size(B);
            FNb = zeros(nJoints, 3, nSamples);
            for s = 1:nSamples
                FNb(:,:,s) = obj.computeFrameNodesOnly(ls, alpha_deg, B(s,:));
            end
        end

        function plotConfigurations(obj,filename, description, genforces, betas, max_segment)
             function plot_model()
                obj.fe.face_alpha = 0.1;

                obj.fe.edge_color='none';
                obj.fe.face_color=[0.2 0.2 0.2]; 
                
                obj.fe.plot(obj.mesh.nodes);
                obj.frameElem.plot(obj.frame_mesh.nodes);
                
                obj.fe.face_alpha = 0.8;
                obj.fe.edge_color=[0.5 0.5 0.5];
                obj.fe.plotSolidSelected(obj.mesh.nodes,obj.const_elems);
                
             end

            obj.fe.face_alpha = 0.1;
            obj.fe.edge_color='none';
            obj.fe.face_color=[0.2 0.2 0.2];  
            ec = obj.fe.edge_color;

            max_elem = max_segment;
            
            figure;
            hold on, axis on; 
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");
            
            tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
 
            % --- Top view ---
            nexttile;
            hold on, axis on; 
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");

            light('Position', [-1 -2 5], 'Style', 'local');
            light('Position', [1 1 5], 'Style', 'infinite');
            plot_model();
            view(2);                % top (xy) view
            axis equal off;
            title('Top view');
        
            % --- Front view ---
            nexttile;
            hold on, axis on; 
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");

            light('Position', [-1 -2 5], 'Style', 'local');
            light('Position', [1 1 5], 'Style', 'infinite');
            plot_model();
            view(0,0);              % front (xz) view
            axis equal off;
            title('Front view');
        
            % --- Side view ---
            nexttile;
            hold on, axis on; 
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");

            light('Position', [-1 -2 5], 'Style', 'local');
            light('Position', [1 1 5], 'Style', 'infinite');
            plot_model();
            view(90,0);             % side (yz) view
            axis equal off;
            title('Side view');
        
            % --- 3D / perspective view ---
            nexttile;
            hold on, axis on; 
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");

            light('Position', [-1 -2 5], 'Style', 'local');
            light('Position', [1 1 5], 'Style', 'infinite');
            plot_model();
            obj.frameElem.plotLocalCS(obj.frame_mesh.nodes); %, 'gamma_deg', betas, 'offsetLocal',[0 0 0.02], 'scale', 0.2);
            obj.frameElem.plotSelected(obj.frame_mesh.nodes, max_elem);
            view(3);                % default 3D
            axis equal off;
            title('3D view');

            sgtitle({ description, [' Q_1= ' sprintf('%.2f',genforces(1)) ', Q_2= ' sprintf('%.2f',genforces(2)) ', Q_3= ' sprintf('%.2f',genforces(3)) ', Q_4= ' sprintf('%.2f',genforces(4)) ', Q_5= ' sprintf('%.2f',genforces(5)) ', Q_6= ' sprintf('%.2f',genforces(6))]}, 'FontWeight','bold');

            exportgraphics(gcf, filename+".png", 'Resolution', 1200); 
            savefig(gcf, filename+".fig");

            figure,hold on, axis on; 
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");
             % --- 3D / perspective view ---
            nexttile;
            hold on, axis on; 
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");

            light('Position', [-1 -2 5], 'Style', 'local');
            light('Position', [1 1 5], 'Style', 'infinite');
            % obj.fe.plot(obj.mesh.nodes);
            % obj.frameElem.plot(obj.frame_mesh.nodes);
            plot_model();
            obj.frameElem.plotLocalCS(obj.frame_mesh.nodes); %, 'gamma_deg', betas, 'offsetLocal',[0 0 0.02], 'scale', 0.2);
            obj.frameElem.plotSelected(obj.frame_mesh.nodes, max_elem);
            view(3);                % default 3D
            axis equal;

            title({ description, [' Q_1= ' sprintf('%.2f',genforces(1)) ', Q_2= ' sprintf('%.2f',genforces(2)) ', Q_3= ' sprintf('%.2f',genforces(3)) ', Q_4= ' sprintf('%.2f',genforces(4)) ', Q_5= ' sprintf('%.2f',genforces(5)) ', Q_6= ' sprintf('%.2f',genforces(6))]}, 'FontWeight','bold', 'FontSize', 28);

            obj.fe.face_alpha = 1.0;
            obj.fe.edge_color=ec;

            exportgraphics(gcf, filename+"_PAN.png", 'Resolution', 1200); 
            savefig(gcf, filename+"_PAN.fig");
        end

        function x = segmentToArm(obj, xOpt)
            nArms=size(obj.mesh.elems,1)/size(xOpt,1)/2;
            x=xOpt;

            for k=1:nArms-1
                x=[x; flip(xOpt); xOpt ];
            end
            x=[x; flip(xOpt)];
        end

        function saveModelMatrices(obj,filename)
            obj.analysis.saveMatrices(filename+"_solid");
            obj.frame_analysis.saveMatrices(filename+"_frame3D");
        end

        function [q_solid, frame_forces] = frameBasedSolver(obj, x)
            % FRAMEBASEDSOLVER  Solve the 3D solid driven by the frame solution
            % via a two-stage condensation scheme.
            %
            % The method is the inverse of the condensation chain used to
            % derive K_beam from K_full (see textbook stages below).  Starting
            % from the frame solution q_frame we expand back to the full solid
            % displacement field q_solid.
            %
            % ----------------------------------------------------------------
            % Stage 1 – Kinematic condensation  (textbook Stage 2, reversed)
            %
            %   Each solid boundary node p is tied to a master frame node m
            %   via the rigid-body transformation matrix T_im:
            %
            %     q_d(p) = T_im * u_m,   where
            %
            %              [ 1  0  0 |  0   Δz  -Δy ]
            %     T_im  =  [ 0  1  0 | -Δz  0    Δx ]   Δ = p − c_m
            %              [ 0  0  1 |  Δy -Δx   0  ]
            %
            %   Evaluated as:  q_d(p) = u_m + θ_m × (p − c_m)
            %
            %   Only the root (frame node 1) and tip (last frame node) act as
            %   kinematic masters.  Interior joints stay free to prevent
            %   double-clamping and the associated kinematic locking.
            %
            % Stage 2 – Static condensation  (textbook Stage 1, reversed)
            %
            %   Partition K_full of the solid into free (i) and prescribed (d):
            %
            %     K_full = [ K_ii  K_id ]
            %              [ K_di  K_dd ]
            %
            %   With zero external load (P_i = 0) and q_d from Stage 1, solve:
            %
            %     K_ii * q_i = −K_id * q_d
            %
            % ----------------------------------------------------------------
            % INPUTS
            %   x  - element density vector [nElems x 1]; defaults to all ones.
            % OUTPUTS
            %   q_solid      - solid nodal displacements [nNodes x 3]
            %   frame_forces - local frame forces [12 x nFrameElems]

            if nargin < 2 || isempty(x) || (isscalar(x) && x == 1)
                x = ones(obj.analysis.getTotalElemsNumber(), 1);
            end

            % ---- Solve frame model to obtain master DOFs u_m ----------------
            nFE_frame = size(obj.frame_mesh.elems, 1);
            obj.frame_analysis.solveWeighted(ones(nFE_frame, 1));
            q_frame = obj.frame_analysis.qnodal;   % [nFrameNodes x 6]: [ux uy uz θx θy θz]

            dof_tr  = obj.frame_analysis.findDOFsIndices(["ux","uy","uz"]);
            dof_rot = obj.frame_analysis.findDOFsIndices(["fix","fiy","fiz"]);

            % =================================================================
            % Stage 1 – Kinematic condensation
            %   Prescribe dependent solid boundary DOFs q_d from frame master
            %   DOFs u_m using the rigid-body relation:
            %     q_d(p) = u_m + θ_m × (p − c_m)   [= T_im * u_m]
            % =================================================================
            nSN         = size(obj.mesh.nodes, 1);
            q_d         = zeros(nSN, 3);     % prescribed (dependent) displacements
            is_d        = false(nSN, 3);     % flags dependent DOFs

            couplingSections = obj.buildCouplingSections();
            if obj.debugCoupling
                obj.reportCouplingSectionDiagnostics(couplingSections);
            end

            nFN          = size(obj.frameNodes, 1);
            masterNodes  = [1, nFN];         % root and tip are the only masters
            ownerMaster  = zeros(nSN, 1);

            for s = 1:numel(couplingSections)
                k = couplingSections(s).frameNode;
                if ~ismember(k, masterNodes)
                    continue;               % interior joints: leave free (no locking)
                end
                p_ids = couplingSections(s).nodes;
                if isempty(p_ids), continue; end

                conflict = ownerMaster(p_ids) ~= 0 & ownerMaster(p_ids) ~= k;
                if any(conflict)
                    error('Overlapping coupling sections: section %d (%s, frame node %d) conflicts with frame node(s) %s.', ...
                        s, couplingSections(s).label, k, mat2str(unique(ownerMaster(p_ids(conflict)))'));
                end
                ownerMaster(p_ids) = k;

                c_m  = obj.frameNodes(k, :);              % master position
                u_m  = q_frame(k, dof_tr);                % master translations
                th_m = q_frame(k, dof_rot);               % master rotations θ_m
                dp   = obj.mesh.nodes(p_ids, :) - c_m;   % (p − c_m) [M×3]

                q_d(p_ids, :)  = repmat(u_m, size(dp,1), 1) + cross(repmat(th_m, size(dp,1), 1), dp, 2);
                is_d(p_ids, :) = true;
            end

            % =================================================================
            % Stage 2 – Static condensation
            %   Assemble K_full, partition into (i) free and (d) prescribed DOFs,
            %   then solve  K_ii * q_i = −K_id * q_d  with P_i = 0.
            % =================================================================
            is_d_full = reshape((logical(obj.analysis.supports) | is_d)', [], 1);
            q_d_full  = reshape(q_d', [], 1);

            [I, J, ~] = obj.analysis.globalMatrixIndices();
            K_vals    = obj.analysis.globalMatrixAggregationWeighted('computeStifnessMatrix', x);
            % K_full = [K_ii  K_id; K_di  K_dd]  (partitioned by is_d_full)

            solver_sc = LinearEquationsSystem(I, J, is_d_full);
            P_i       = zeros(numel(q_d_full), 1);   % no external load on free DOFs
            q_full    = solver_sc.solvePq(K_vals, P_i, q_d_full);

            % Restore prescribed values q_d into the full solution vector
            q_full(solver_sc.supdofs) = q_d_full(solver_sc.supdofs);
            if obj.debugCoupling
                fprintf('Stage 2 restore error: %.3e\n', ...
                    max(abs(q_full(solver_sc.supdofs) - q_d_full(solver_sc.supdofs))));
            end
            assert(max(abs(q_full(solver_sc.supdofs) - q_d_full(solver_sc.supdofs))) < 1.0e-12, ...
                'Prescribed DOFs were not restored correctly after static condensation.');

            % ---- Store results and compute element stresses -----------------
            obj.analysis.qfem   = q_full;
            obj.analysis.qnodal = obj.analysis.fromFEMVector(q_full);
            obj.qnodal_solid    = obj.analysis.qnodal;
            obj.analysis.computeElementResults(x);

            q_solid = obj.qnodal_solid;
            [~, frame_forces] = obj.frameElem.computeResults(obj.frame_mesh.nodes, q_frame);
        end


        function couplingSections = buildCouplingSections(obj)
            % The exact section-node IDs are captured during segment
            % generation and mapped to global IDs after each merge. This
            % avoids reconstructing coupling sets from exterior shell nodes.
            couplingSections = obj.couplingSections;
        end

        function edgeSets = buildEdgeNodeSets(obj)
            % Legacy compatibility wrapper: return one node set per frame
            % node by grouping the exact per-section coupling sets. The
            % solver itself uses buildCouplingSections() to avoid re-merging
            % different physical sections into one Dirichlet set.
            nFN = size(obj.frameNodes, 1);
            edgeSets = cell(nFN, 1);
            couplingSections = obj.buildCouplingSections();
            for k = 1:nFN
                edgeSets{k} = zeros(0, 1, 'int32');
            end
            for s = 1:numel(couplingSections)
                k = couplingSections(s).frameNode;
                edgeSets{k} = unique([edgeSets{k}; couplingSections(s).nodes], 'stable');
            end
        end

        function plotCouplingSections(obj, sectionIds)
            if nargin < 2 || isempty(sectionIds)
                sectionIds = 1:numel(obj.couplingSections);
            end
            figure;
            hold on;
            axis on;
            daspect([1 1 1]);
            xlabel("x");
            ylabel("y");
            zlabel("z");
            obj.fe.face_alpha = 0.1;
            obj.fe.edge_color = 'none';
            obj.fe.plot(obj.mesh.nodes);
            cmap = lines(max(numel(sectionIds), 1));
            for i = 1:numel(sectionIds)
                s = sectionIds(i);
                inodes = obj.couplingSections(s).nodes;
                scatter3(obj.mesh.nodes(inodes,1), obj.mesh.nodes(inodes,2), obj.mesh.nodes(inodes,3), ...
                    40, repmat(cmap(i,:), numel(inodes), 1), 'filled');
            end
        end

        function gids = mapNodesByCoordinates(obj, allNodes, nodeCoords, tol)
            roundedAll = round(allNodes * tol);
            roundedCoords = round(nodeCoords * tol);
            [tf, gids] = ismember(roundedCoords, roundedAll, 'rows');
            assert(all(tf), 'Failed to map section nodes after mesh merge.');
            gids = int32(gids(:));
        end

        function registerCouplingSection(obj, frameNode, adjacentFrameNode, nodeIds, label)
            sec.frameNode = frameNode;
            sec.adjacentFrameNode = adjacentFrameNode;
            sec.nodes = unique(int32(nodeIds(:)), 'stable');
            sec.label = char(label);
            obj.couplingSections(end+1) = sec;
        end

        function reportCouplingSectionDiagnostics(obj, couplingSections)
            fprintf('Coupling section diagnostics:\n');
            for s = 1:numel(couplingSections)
                sec = couplingSections(s);
                ck = obj.frameNodes(sec.frameNode, :);
                normal = obj.frameNodes(sec.adjacentFrameNode, :) - ck;
                normal = normal / norm(normal);
                dp = obj.mesh.nodes(sec.nodes, :) - ck;
                axial = dp * normal';
                nComp = obj.countCouplingSectionComponents(sec.nodes);
                fprintf('  section %d (%s): frame node %d, adjacent %d, nodes=%d, axial span=[%.3e, %.3e], components=%d\n', ...
                    s, sec.label, sec.frameNode, sec.adjacentFrameNode, numel(sec.nodes), ...
                    min(axial), max(axial), nComp);
            end
        end

        function nComp = countCouplingSectionComponents(obj, nodeIds)
            allEdges = obj.fe.multiObjectList(obj.fe.shapeFn.edges);
            edgeMask = all(ismember(allEdges, nodeIds), 2);
            sectionEdges = unique(sort(allEdges(edgeMask, :), 2), 'rows');
            if isempty(sectionEdges)
                nComp = numel(nodeIds);
                return;
            end

            sectionNodeIds = unique(sectionEdges(:));
            nComp = 0;
            visited = false(numel(sectionNodeIds), 1);
            adjacency = cell(numel(sectionNodeIds), 1);
            [~, locA] = ismember(sectionEdges(:,1), sectionNodeIds);
            [~, locB] = ismember(sectionEdges(:,2), sectionNodeIds);
            for e = 1:size(sectionEdges, 1)
                adjacency{locA(e)} = [adjacency{locA(e)}; locB(e)];
                adjacency{locB(e)} = [adjacency{locB(e)}; locA(e)];
            end

            for n = 1:numel(sectionNodeIds)
                if visited(n)
                    continue;
                end
                nComp = nComp + 1;
                queue = n;
                visited(n) = true;
                while ~isempty(queue)
                    current = queue(1);
                    queue(1) = [];
                    neighbors = adjacency{current};
                    for m = neighbors(:)'
                        if ~visited(m)
                            visited(m) = true;
                            queue(end+1) = m; %#ok<AGROW>
                        end
                    end
                end
            end

            isolatedNodes = setdiff(nodeIds(:), sectionNodeIds(:));
            nComp = nComp + numel(isolatedNodes);
        end

        function applyJointFillets(obj, filletR, xEnds, R_outer)
            % applyJointFillets  Round the outer-wall ridge at each interior bevel joint.
            %
            % For each interior joint (not root, not tip), nodes within the
            % fillet zone — within filletR of both the bevel plane and the
            % outer surface — are displaced inward radially to follow a
            % circular arc of radius filletR.  Nodes at the bevel plane
            % (d_plane=0) are not moved; the arc dips inward symmetrically
            % on both sides, replacing the sharp ridge with a smooth convex
            % fillet.  Only nodes that currently sit outside the arc
            % (r_node > r_target) are touched, so inner-wall nodes and
            % mid-body nodes are unaffected.
            if filletR <= 0
                return;
            end
            Th = R_outer - obj.r;
            if filletR > 0.5 * Th
                warning('ManipulatorModel3D:filletR', ...
                    'filletR (%.3g m) > 0.5*wallThickness (%.3g m); clamping to avoid element inversion.', ...
                    filletR, 0.5 * Th);
                filletR = 0.5 * Th;
            end

            nFrameNodes = size(xEnds, 1);
            nodes = obj.mesh.nodes;

            for k = 2:nFrameNodes - 1
                p_joint = xEnds(k, :);

                d_in  = p_joint - xEnds(k-1, :);
                d_out = xEnds(k+1, :) - p_joint;
                if norm(d_in) < 1e-12 || norm(d_out) < 1e-12, continue; end
                d_in  = d_in  / norm(d_in);
                d_out = d_out / norm(d_out);

                % Bevel plane normal: bisector of the two adjacent tube axes.
                n_bev = d_in + d_out;
                if norm(n_bev) < 1e-10, continue; end
                n_bev = n_bev / norm(n_bev);

                % Signed distance of every node from the bevel plane.
                dp      = nodes - p_joint;                       % [N x 3]
                d_plane = dp * n_bev';                           % [N x 1]

                % Radial distance from the line through p_joint along n_bev.
                axial_comp = d_plane .* n_bev;                   % [N x 3]
                radial_vec = dp - axial_comp;                    % [N x 3]
                r_node     = sqrt(sum(radial_vec .^ 2, 2));      % [N x 1]

                % Fillet zone: near bevel plane AND near outer surface.
                in_zone = abs(d_plane) < filletR & r_node > R_outer - filletR;
                if ~any(in_zone), continue; end

                % Arc: r_target = R_outer - filletR + sqrt(filletR^2 - d_plane^2).
                % Move only nodes that currently protrude outside the arc.
                idx   = find(in_zone);
                d_sq  = d_plane(idx) .^ 2;
                r_tgt = R_outer - filletR + sqrt(max(0, filletR^2 - d_sq));

                move  = r_node(idx) > r_tgt + 1e-14;
                if ~any(move), continue; end

                imove = idx(move);
                scale = r_tgt(move) ./ r_node(imove);
                nodes(imove, :) = p_joint + axial_comp(imove, :) + ...
                    radial_vec(imove, :) .* scale;
            end

            obj.mesh.nodes = nodes;
        end

    end
end

classdef ManipulatorModel3D < handle
 
    
    properties
        analysis, frame_analysis, fixedEdgeSelector, alpha, mesh, frame_mesh, elems, fe, xEnd, frameNodes, frameElem
        const_elems, upper_nodes, loadSurfaceNodes, fixedSurfaceNodes, halfSegmentNelems, use_offset;
    end
    
    methods                       
        function obj = ManipulatorModel3D(E,nu,ls,R,r, res, res_th, alpha, betas, ShapeFn, use_offset, Pz)
            alpha=alpha*pi/180;
            betas=betas*pi/180;

            if use_offset
                obj.use_offset = 1;
            else
                obj.use_offset = 0;
            end

            obj.mesh=Mesh();
            obj.geterateManipulator( ls, R, r, res, res_th, alpha, betas, ShapeFn );
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
            
            obj.analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 Pz/A] ));
            obj.analysis.fixNodes( fixedEdgeSelector, [ "uz"] );
            obj.analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);

            frame_elems = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
            obj.frame_mesh=Mesh();
            obj.frame_mesh.nodes=obj.frameNodes;
            obj.frame_mesh.elems = frame_elems;
            
            obj.frameElem=Frame3D(frame_elems,E,nu,R,r);
            obj.frame_analysis = LinearElasticityWeighted( obj.frameElem, obj.frame_mesh, false );
            obj.frame_analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 0 0 0 0]);
            obj.frame_analysis.loadClosestNode(obj.frame_mesh.nodes(end,:), ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 -Pz 0 0 0] );

        end

        function ns = determineSegment( obj, elem )
            nhs = round(elem ./ obj.halfSegmentNelems );  
            ns = ceil((nhs - 1) / 2) + 1;
        end

        function geterateManipulator(obj, ls, R, r, res, res_th, alpha, betas, sf )
            Th = R - r;
        
            % --- thickness resolution (radial) ---
            resTh = max(1, round(res_th));
        
            % --- choose circum + length resolution consistently with resTh ---
            % (your original logic: scale with Th so elements stay ~ isotropic)
            resCirc = max(3, round(2*pi*R/Th * resTh));
            resLen  = max(1, round(ls/Th * resTh));
        
            % rotations...
            c=cos(alpha); s=sin(alpha);
            rotCut=[c 0 s; 0 1 0; -s 0 c]';
            c=cos(2*alpha); s=sin(2*alpha);
            rotCut2=[c 0 s; 0 1 0; -s 0 c]';
            c=cos(betas(1)); s=sin(betas(1));
            rotBeta=[c -s 0; s c 0; 0 0 1]';
        
            phase = -betas(1)*obj.use_offset;
        
            % --- first segment ---
            [obj.mesh, obj.elems] = obj.generateSegment2a(R, r, ls, phase, rotCut, sf, resTh, resCirc, resLen);
        
            obj.halfSegmentNelems = size(obj.elems,1);
        
            % elements in ONE end slice (k = last along length):
            sliceElems = resTh * resCirc;
            ne = size(obj.elems,1);
            obj.const_elems = (ne - sliceElems + 1 : ne).';
        
            obj.mesh.nodes = obj.mesh.nodes * rotBeta;
            selector = Selector(@(x)(x(:,3) < 1.0E-4));
            obj.fixedSurfaceNodes = selector.select(obj.mesh.nodes);
        
            prevRot = rotCut * rotBeta;
            obj.xEnd = [0 0 ls];
            xEnds = [[0 0 0]; obj.xEnd];
        
            for k = 2:length(betas)
                c=cos(betas(k)); s=sin(betas(k));
                rotBeta=[c -s 0; s c 0; 0 0 1]';
        
                phase = obj.use_offset * (phase - betas(k));
        
                [mesh1, elems1] = obj.generateSegment2b(R, r, ls, phase, rotCut, sf, resTh, resCirc, resLen);
                [mesh2, elems2] = obj.generateSegment2a(R, r, ls, phase, rotCut, sf, resTh, resCirc, resLen);
        
                if k < length(betas)
                    elems1 = [elems2; mesh1.merge(mesh2.nodes, elems2)];
                end
        
                mesh1.nodes = (mesh1.nodes + [0 0 ls]) * rotCut * rotBeta * prevRot + obj.xEnd;
                obj.elems = [obj.elems; obj.mesh.merge(mesh1.nodes, elems1)];
        
                % update const_elems for this newly appended part:
                ne = size(obj.elems,1);
                obj.const_elems = [obj.const_elems; (ne - sliceElems + 1 : ne).'];
        
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
        end


        function [mesh, elems] = generateSegment2a(obj, R, r, ls, phase, Redge, sf, resTh, resCirc, resLen)
            Th = R - r;
        
            mesh = Mesh();
            elems = mesh.addRectMesh3D( ...
                R - Th, phase, 0, ...
                Th, 2*pi, 1, ...
                resTh, resCirc, resLen, sf.localNodes);
        
            elems = mesh.transformToCylindrical3D([0 0]);
        
            nodesR1 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)] * Redge;
            nodes   = [nodesR1(:,1), nodesR1(:,2), mesh.nodes(:,3) .* (nodesR1(:,3) + ls)];
            mesh.nodes = nodes;
            elems = mesh.elems;
        end


        function [mesh, elems] = generateSegment2b(obj, R, r, ls, phase, Redge, sf, resTh, resCirc, resLen)
            Th = R - r;
        
            mesh = Mesh();
            elems = mesh.addRectMesh3D( ...
                R, phase, 0, ...
                -Th, -2*pi, 1, ...
                resTh, resCirc, resLen, sf.localNodes);
        
            elems = mesh.transformToCylindrical3D([0 0]);
        
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

            title({ description, [' Q_1= ' sprintf('%.2f',genforces(1)) ', Q_2= ' sprintf('%.2f',genforces(2)) ', Q_3= ' sprintf('%.2f',genforces(3)) ', Q_4= ' sprintf('%.2f',genforces(4)) ', Q_5= ' sprintf('%.2f',genforces(5)) ', Q_6= ' sprintf('%.2f',genforces(6))]}, 'FontWeight','bold');

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
        
    end
end


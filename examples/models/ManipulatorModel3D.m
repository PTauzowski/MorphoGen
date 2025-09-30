classdef ManipulatorModel3D < handle
 
    
    properties
        analysis, frame_analysis, fixedEdgeSelector, alpha, mesh, frame_mesh, elems, fe, xEnd, frameNodes, frameElem
        const_elems, upper_nodes, loadSurfaceNodes, fixedSurfaceNodes, halfSegmentNelems, use_offset;
    end
    
    methods                       
        function obj = ManipulatorModel3D(E,nu,ls,R,r,res, alpha, betas, ShapeFn, use_offset)
            alpha=alpha*pi/180;
            betas=betas*pi/180;

            if use_offset
                obj.use_offset = 1;
            else
                obj.use_offset = 0;
            end

            obj.mesh=Mesh();
            obj.geterateManipulator( ls, R, r, res, alpha, betas, ShapeFn );
            obj.fe = SolidElasticElem( ShapeFn, obj.elems );

            obj.fe.props.h=1;
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material);

            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
           
            fixedEdgeSelector = Selector( obj.fixedSurfaceNodes );
            loadedFaceSelector = Selector( obj.loadSurfaceNodes );
            
            obj.analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 -10] ));
            obj.analysis.fixNodes( fixedEdgeSelector, [ "uz"] );
            obj.analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);

            frame_elems = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8];
            obj.frame_mesh=Mesh();
            obj.frame_mesh.nodes=obj.frameNodes;
            obj.frame_mesh.elems = frame_elems;
            
            obj.frameElem=Frame3D(frame_elems,E,0.02,0.8*E,0.0004,0.0004,0.003);
            obj.frame_analysis = LinearElasticityWeighted( obj.frameElem, obj.frame_mesh, false );
            obj.frame_analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 0 0 0 0]);
            obj.frame_analysis.loadClosestNode(obj.frame_mesh.nodes(end,:), ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 -1 0 0 0] );

        end

        function ns = determineSegment( obj, elem )
            nhs = round(elem ./ obj.halfSegmentNelems );  
            ns = ceil((nhs - 1) / 2) + 1;
        end

        function geterateManipulator(obj, ls, R, r, res, alpha, betas, sf )
            %obj.elems = obj.mesh.merge(mesh.nodes,obj.elems)
           Th=R-r;
           resTh=1;
           resCirc=round(2*pi*R/Th*resTh);
           c=cos(alpha);
           s=sin(alpha);
           rotCut=[c 0 s; 0 1 0; -s 0 c]';
           c=cos(2*alpha);
           s=sin(2*alpha);
           rotCut2=[c 0 s; 0 1 0; -s 0 c]';
           c=cos(betas(1));
           s=sin(betas(1));
           rotBeta=[c -s 0; s c 0; 0 0 1]';
            
            phase=-betas(1)*obj.use_offset;
            [obj.mesh, obj.elems]=obj.generateSegment2a(R, r, ls, res, phase, rotCut, sf);
            obj.halfSegmentNelems=size(obj.elems,1);
            obj.const_elems=(size(obj.elems,1)-1*resCirc-1:size(obj.elems,1))';
            last_const_elems=obj.elems(size(obj.elems,1)-1*resCirc-1:size(obj.elems,1),:);
            obj.mesh.nodes=obj.mesh.nodes*rotBeta;
            selector = Selector( @(x)( (x(:,3) < 1.0E-4) ) );
            obj.fixedSurfaceNodes = selector.select( obj.mesh.nodes );
            prevRot=rotCut*rotBeta;
            obj.xEnd=[0 0 ls];
            xEnds=[ [0 0 0]; obj.xEnd];
           
            for k=2:length(betas)
                c=cos(betas(k));
                s=sin(betas(k));
                rotBeta=[c -s 0; s c 0; 0 0 1]';
                phase = obj.use_offset * ( phase-betas(k) );
                [mesh1, elems1]=obj.generateSegment2b(R, r, ls, res, phase, rotCut, sf);
                [mesh2, elems2]=obj.generateSegment2a(R, r, ls, res, phase, rotCut, sf);
                if k<length(betas)
                    elems1=[elems2; mesh1.merge(mesh2.nodes,elems2)];                    
                end
                
                mesh1.nodes=(mesh1.nodes+[0 0 ls])*rotCut*rotBeta*prevRot+obj.xEnd;  
                obj.elems =[ obj.elems; obj.mesh.merge(mesh1.nodes, elems1) ];
                % last_const_elems=obj.elems(size(obj.elems,1)-5*res-1:size(obj.elems,1),:);
                % obj.const_elems=[ obj.const_elems; (size(obj.elems,1)-5*res:size(obj.elems,1))' ];
                last_const_elems=obj.elems(size(obj.elems,1)-1*resCirc-1:size(obj.elems,1),:);
                obj.const_elems=[ obj.const_elems; (size(obj.elems,1)-1*resCirc:size(obj.elems,1))' ];
                if k<length(betas)
                    obj.xEnd=obj.xEnd+[0 0 2*ls]*rotCut*rotBeta*prevRot;
                else
                    obj.xEnd=obj.xEnd+[0 0 ls]*rotCut*rotBeta*prevRot;
                end
                xEnds=[ xEnds; obj.xEnd ];
                prevRot=rotCut2*rotBeta*prevRot;
            end
            tNodes = (obj.mesh.nodes-repmat(obj.xEnd,size(obj.mesh.nodes,1),1))*prevRot'*rotCut;
            sNodes = abs(tNodes(:,3))<1.0E-4;
            obj.loadSurfaceNodes = sNodes;
            obj.frameNodes=xEnds;
            %plot3(obj.mesh.nodes(sNodes,1),obj.mesh.nodes(sNodes,2),obj.mesh.nodes(sNodes,3),LineStyle="none",Marker="*",Color='r');
            %line(tNodes(:,1),tNodes(:,2),tNodes(:,3),Marker="o",Color='r');
        end

        function [mesh, elems] = generateSegment2a(obj, R, r, ls, res, phase, Redge, sf)
            Th=R-r;
            % resLen=res;
            % resCirc = ceil(resLen/ls*2*pi*R);
            % resTh = ceil(resLen/ls*Th);

            resTh=1;
            resCirc=round(2*pi*R/Th*resTh);
            resLen=round(ls/Th*resTh);
            
            mesh = Mesh();
            elems = mesh.addRectMesh3D( R-Th, phase, 0, Th, 2*pi, 1, resTh, resCirc, resLen, sf.localNodes );      
            elems = mesh.transformToCylindrical3D( [0 0] );
           
            nodesR1 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)]*Redge;
            %nodesR2 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)]*Rot2;
            nodes= [nodesR1(:,1) nodesR1(:,2) mesh.nodes(:,3) .* (nodesR1(:,3)+ls)];
            mesh.nodes=nodes;   
            elems=mesh.elems;
        end

        function [mesh, elems] = generateSegment2b(obj, R, r, ls, res, phase, Redge, sf)
            Th=R-r;
            % resLen=res;
            % resCirc = ceil(resLen/ls*2*pi*R);
            % resTh = ceil(resLen/ls*Th);

            resTh=1;
            resCirc=round(2*pi*R/Th*resTh);
            resLen=round(ls/Th*resTh);
            
            mesh = Mesh();
            elems = mesh.addRectMesh3D( R, phase, 0, -Th, -2*pi, 1, resTh, resCirc, resLen, sf.localNodes );      
            elems = mesh.transformToCylindrical3D( [0 0] );
           
            nodesR1 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)]*Redge';
            %nodesR2 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)]*Rot2;
            nodes= [nodesR1(:,1) nodesR1(:,2) (1-mesh.nodes(:,3)) .* (nodesR1(:,3)-ls)];
            mesh.nodes=nodes; 
            elems=mesh.elems;
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
            figure
            hold on, axis on; 
            daspect([1 1 1]);
            obj.analysis.felems{1}.plot(obj.mesh.nodes);
            obj.analysis.plotCurrentLoad();
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

        function plotConfigurations(obj,filename, description, genforces, max_segment)
            % figure
            % hold on, axis on; 
            % daspect([1 1 1]);
            % xlabel("x");
            % ylabel("y");
            % zlabel("z");
            % 
            % light('Position', [-1 -2 5], 'Style', 'local');
            % light('Position', [1 1 5], 'Style', 'infinite');
            % ax = gca; 
            % 

            max_elem = max_segment;

            obj.fe.face_alpha = 0.2;
            ec = obj.fe.edge_color;
            obj.fe.edge_color='none';
            obj.fe.face_color=[0.2 0.2 0.2];
            
            figure;
            
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
            obj.fe.plot(obj.mesh.nodes);
            obj.frameElem.plot(obj.frame_mesh.nodes);
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
            obj.fe.plot(obj.mesh.nodes);
            obj.frameElem.plot(obj.frame_mesh.nodes);
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
            obj.fe.plot(obj.mesh.nodes);
            obj.frameElem.plot(obj.frame_mesh.nodes);
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
            obj.fe.plot(obj.mesh.nodes);
            obj.frameElem.plot(obj.frame_mesh.nodes);
            obj.frameElem.plotLocalCS(obj.frame_mesh.nodes, 0.2, 0.02);
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
            obj.fe.plot(obj.mesh.nodes);
            obj.frameElem.plot(obj.frame_mesh.nodes);
            obj.frameElem.plotLocalCS(obj.frame_mesh.nodes, 0.2,0.02);
            obj.frameElem.plotSelected(obj.frame_mesh.nodes, max_elem);
            view(3);                % default 3D
            axis equal;

            title({ description, [' Q_1= ' sprintf('%.2f',genforces(1)) ', Q_2= ' sprintf('%.2f',genforces(2)) ', Q_3= ' sprintf('%.2f',genforces(3)) ', Q_4= ' sprintf('%.2f',genforces(4)) ', Q_5= ' sprintf('%.2f',genforces(5)) ', Q_6= ' sprintf('%.2f',genforces(6))]}, 'FontWeight','bold');

            obj.fe.face_alpha = 1.0;
            obj.fe.edge_color=ec;

            exportgraphics(gcf, filename+"_PAN.png", 'Resolution', 1200); 
            savefig(gcf, filename+"_PAN.fig");
        end
        
    end
end


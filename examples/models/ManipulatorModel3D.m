classdef ManipulatorModel3D < handle
 
    
    properties
        analysis, fixedEdgeSelector, alpha, mesh, elems, fe, xEnd, const_elems, upper_nodes, loadSurfaceNodes, fixedSurfaceNodes, frameNodes;
    end
    
    methods                       
        function obj = ManipulatorModel3D(E,nu,ls,R,r,res, alpha, betas, ShapeFn)
            alpha=alpha*pi/180;
            betas=betas*pi/180;

            obj.mesh=Mesh();
            obj.geterateManipulator2( ls, R, r, res, alpha, betas, ShapeFn );
            obj.fe = SolidElasticElem( ShapeFn, obj.elems );

            obj.fe.props.h=1;
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material);

            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
            %analysis = LinearElasticityWeighted( model.fe, model.mesh, false );
            %problem = LinearElasticity( fe, mesh );
            fixedEdgeSelector = Selector( obj.fixedSurfaceNodes );
            loadedFaceSelector = Selector( obj.loadSurfaceNodes );
            
            obj.analysis.elementLoadSurfaceIntegral( "global", loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [0 0 -10] ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
            obj.analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0]);

        end

        function geterateManipulator(obj, ls, R, r, res, alpha, betas, sf )
            th=R-r;
            %obj.elems = obj.mesh.merge(mesh.nodes,obj.elems);
            x0=[0 0 0];
            x1=[0 0 ls/2];
            c=cos(alpha);
            s=sin(alpha);
            Ryp=[c 0 s; 0 1 0; -s 0 c];
            c=cos(-alpha);
            s=sin(-alpha);
            Rym=[c 0 s; 0 1 0; -s 0 c];
            Rot1=eye(3);
            Rot2=Ryp;
            Rotk=Rot2;
            xs=x1;
            [obj.mesh, obj.elems]=obj.generateSegment(R, r, ls, res, x0, x1, 0, Rot1, Rot2, sf);
            obj.const_elems=(size(obj.elems,1)-2*res-1:size(obj.elems,1))';
            phase = 0;
            for k=1:length(betas)
                c=cos(betas(k));
                s=sin(betas(k));
                Rb=[c -s 0; s c 0; 0 0 1];
                Rotm=Rym*Rb*Rotk;
                phase=phase+betas(k);
                xm=[0 0 ls/2]*Rotm+xs;
                [mesh, elems]=obj.generateSegment(R, r, ls, res, xs, xm, phase, Rb*Rotk, Rotm, sf);
                obj.elems = [obj.elems; obj.mesh.merge(mesh.nodes, elems)];
                obj.const_elems=[ obj.const_elems; (size(obj.elems,1)-2*res:size(obj.elems,1))' ];
                xs=[0 0 ls]*Rotm+xs;
                Rotk=Rym*Rotm;
                [mesh, elems]=obj.generateSegment(R, r, ls, res, xm, xs, phase, Rotm, Rotk, sf);
                obj.elems = [obj.elems; obj.mesh.merge(mesh.nodes, elems)];
                last_const_elems=obj.elems(size(obj.elems,1)-2*res-1:size(obj.elems,1),:);
                obj.const_elems=[ obj.const_elems; (size(obj.elems,1)-2*res:size(obj.elems,1))' ];
            end
            ne=size(last_const_elems,1);
            nodesCount=zeros(size(obj.mesh.nodes,1),1);
            for k=1:ne
                nodesCount(last_const_elems(k,:))=nodesCount(last_const_elems(k,:))+1;
            end
            obj.loadSurfaceNodes=nodesCount==2;
        end

        function [mesh, elems] = generateSegment(obj, R, r, ls, res, x1, x2, phase, Rot1, Rot2, sf)
            Th=R-r;
            resLen=res;
            resCirc = ceil(resLen/ls*2*pi*R);
            resTh = ceil(resLen/ls*Th);
            
            mesh = Mesh();
            elems = mesh.addRectMesh3D( R-Th, phase, 0, Th, 2*pi, 1, resTh, resCirc, resLen, sf.localNodes );      
            elems = mesh.transformToCylindrical3D( [0 0], elems );
           
            nodesR1 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)]*Rot1+x1;
            nodesR2 = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)]*Rot2+x2;
            nodes= (1-mesh.nodes(:,3)) .* nodesR1 + mesh.nodes(:,3) .* nodesR2;
            mesh.nodes=nodes;            
        end

        function geterateManipulator2(obj, ls, R, r, res, alpha, betas, sf )
            %obj.elems = obj.mesh.merge(mesh.nodes,obj.elems)
            
           c=cos(alpha);
           s=sin(alpha);
           rotCut=[c 0 s; 0 1 0; -s 0 c]';
           c=cos(2*alpha);
           s=sin(2*alpha);
           rotCut2=[c 0 s; 0 1 0; -s 0 c]';
           c=cos(betas(1));
           s=sin(betas(1));
           rotBeta=[c -s 0; s c 0; 0 0 1]';
            
            phase=-betas(1);
            [obj.mesh, obj.elems]=obj.generateSegment2a(R, r, ls, res, phase, rotCut, sf);
            obj.const_elems=(size(obj.elems,1)-2*res-1:size(obj.elems,1))';
            last_const_elems=obj.elems(size(obj.elems,1)-2*res-1:size(obj.elems,1),:);
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
                phase = phase-betas(k);
                [mesh1, elems1]=obj.generateSegment2b(R, r, ls, res, phase, rotCut, sf);
                [mesh2, elems2]=obj.generateSegment2a(R, r, ls, res, phase, rotCut, sf);
                if k<length(betas)
                    elems1=[elems2; mesh1.merge(mesh2.nodes,elems2)];                    
                end
                
                mesh1.nodes=(mesh1.nodes+[0 0 ls])*rotCut*rotBeta*prevRot+obj.xEnd;  
                obj.elems =[ obj.elems; obj.mesh.merge(mesh1.nodes, elems1) ];
                last_const_elems=obj.elems(size(obj.elems,1)-2*res-1:size(obj.elems,1),:);
                obj.const_elems=[ obj.const_elems; (size(obj.elems,1)-2*res:size(obj.elems,1))' ];
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
            elems = mesh.addRectMesh3D( R-Th, phase, 0, Th, 2*pi, 1, resTh, resCirc, resLen, sf.localNodes );      
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
        
    end
end


classdef ManipulatorModel3D < FEModel
 
    
    properties
        analysis, fixedEdgeSelector, alpha, elems, fe, const_elems, upper_nodes,loadSurfaceNodes;
    end
    
    methods                       
        function obj = ManipulatorModel3D(ls,R,r,res, alpha, betas, ShapeFn)
            alpha=alpha*pi/180;
            betas=betas*pi/180;

            obj.mesh=Mesh();
            obj.geterateManipulator( ls, R, r, res, alpha, betas, ShapeFn );
            obj.fe = SolidElasticElem( ShapeFn, obj.elems );

            obj.fe.props.h=1;
            material = SolidMaterial('mat1');
            material.setElasticIzo(210.0E9, 0.3);
            material.setElasticIzoGrad();

            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
            %problem = LinearElasticity( fe, mesh );
            fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.001 ) );
            loadedFaceSelector = Selector( @(x)( abs(x(:,3)- Length) < 0.001 ) );
            constElemsSelector =  Selector( @(x)( (x(:,3) < 0.05 * Length ) ) & (x(:,3) > 0.96 * Length ) );
            %obj.setX(ones(obj.analysis.getTotalElemsNumber(),1));
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

        function xnew=rotatePoints(x0,R,x)
            xnew=(x-x0)*R+x0;
        end       
        
    end
end


classdef ManipulatorModel3D < FEModel
 
    
    properties
        analysis, fixedEdgeSelector, alpha, elems, fe;
    end
    
    methods                       
        function obj = ManipulatorModel3D(ls,R,r,res, alpha, betas, ShapeFn)
            alpha=alpha*pi/180;
            betas=betas*pi/180;
            Th = R-r;
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
            [mesh, obj.elems]=obj.generateSegment(ls, R, r, res, alpha, sf);
            %obj.elems = obj.mesh.merge(mesh.nodes,obj.elems);
            x0=[0 0 ls];
            c=cos(alpha);
            s=sin(alpha);
            Ry=[c 0 s; 0 1 0; -s 0 c];
            Rot=R;
            for k=1:length(betas)
                [mesh1, elems1]=obj.generateSegment(ls, R, r, res, alpha, sf);
                [mesh2, elems2]=obj.generateSegment(ls, R, r, res, alpha, sf);
                c=cos(pi);
                s=sin(pi);
                Ry=[c 0 s; 0 1 0; -s 0 c];
                Rz=[c -s 0; s c 0; 0 0 1];
                mesh1.nodes=mesh1.nodes*Ry;
                mesh1.nodes=mesh1.nodes*Rz;
                elems=[ elems1; mesh1.merge(mesh2.nodes, elems2)];
                mesh1.nodes=mesh1.nodes*Ry;
                mesh1.nodes = mesh1.nodes+[0 0 ls];
                c=cos(-alpha);
                s=sin(-alpha);
                Ry=[c 0 s; 0 1 0; -s 0 c];
                c=cos(betas(k));
                s=sin(betas(k));
                 Rz=[c -s 0; s c 0; 0 0 1];
                mesh1.nodes=mesh1.nodes*Rz;
                obj.elems=[ obj.elems; obj.mesh.merge(mesh1.nodes, elems) ];
                x0=x0+[0 0 2*ls];
            end
        end

        function [mesh, elems] = generateSegment(obj,ls, R, r, res, alpha, sf)
            Th=R-r;
            resLen=res;
            resCirc = ceil(resLen/ls*2*pi*R);
            resTh = ceil(resLen/ls*Th);
            mesh = Mesh();
            elems = mesh.addRectMesh3D( R-Th, 0, 0, Th, 2*pi, ls, resTh, resCirc, resLen, sf.localNodes );
            elems = mesh.transformToCylindrical3D( [0 0], elems );
            c=cos(alpha);
            s=sin(alpha);
            Ry=[c 0 s; 0 1 0; -s 0 c];
            nodesR = [mesh.nodes(:,1) mesh.nodes(:,2) 0*mesh.nodes(:,3)]*Ry;
            nodes  = [ (1-mesh.nodes(:,3)./ls) .* mesh.nodes(:,1) + mesh.nodes(:,3)./ls .* nodesR(:,1)...
                       (1-mesh.nodes(:,3)./ls) .* mesh.nodes(:,2) + mesh.nodes(:,3)./ls .* nodesR(:,2)... 
                         mesh.nodes(:,3)./ls .* (ls+nodesR(:,3))];
            mesh.nodes=nodes;

        end
        
        
    end
end


classdef ColumnModel3D < ModelLinear
    
    methods
        function obj = ColumnModel3D(sf,b,h,l,nl,E,nu,rho,P,xp)
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addRectMesh3D( 0, 0, 0, b, h, l, round(nl*b/l), round(nl*h/l), nl, sf.localNodes );
            obj.x=ones(size(obj.mesh.elems,1),1);
            fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.001 ) );
            loadedEdgeSelector = Selector( @(x)( abs(x(:,3)-l) < 0.001 ) );

            obj.fe=SolidElasticElem( sf, obj.mesh.elems );
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.rho=rho;
            obj.fe.setMaterial( material );            
            %obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true ); 
            obj.analysis = SecondOrderElasticityWeighted( obj.fe, obj.mesh, false );
            obj.analysis.loadClosestNode(xp, ["ux" "uy" "uz"], P );
            obj.analysis.elementLoadSurfaceIntegral( "global", loadedEdgeSelector,  ["ux" "uy" "uz"], @(x)( x*0 + P ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
            % obj.analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0] );
            % obj.analysis.fixNodes( loadedEdgeSelector, ["ux" "uy"] );
            % obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
            %obj.analysis.fixClosestNode( [b/2 h/2 0], ["ux" "uy" "uz"], [0 0 0] );
            obj.analysis.printProblemInfo();
            obj.plotModel();
        end
    end
end



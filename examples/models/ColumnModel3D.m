classdef ColumnModel3D < ModelLinear
    
    methods
        function obj = ColumnModel3D(sf,l,h,nh,E,nu,P,xp)
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addRectMesh3D( 0, 0, 0, l, l, h, round(nh*l/h), round(nh*l/h), nh, sf.localNodes );
            obj.x=ones(size(obj.mesh.elems,1),1);
            fixedEdgeSelector = Selector( @(x)( abs(x(:,3)) < 0.001 ) );
            loadedEdgeSelector = Selector( @(x)( abs(x(:,3)-h) < 0.001 ) );

            obj.fe=SolidElasticElem( sf, obj.mesh.elems );
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );            
            %obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true ); 
            obj.analysis = SecondOrderElasticityWeighted( obj.fe, obj.mesh, false );
            %obj.analysis.loadClosestNode(xp, ["ux" "uy" "uz"], P );
            obj.analysis.elementLoadSurfaceIntegral( "global", loadedEdgeSelector,  ["ux" "uy" "uz"], @(x)( x*0 + P ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy" "uz"] );
            %obj.analysis.fixClosestNode( [0 0 0], ["ux" "uy" "uz"], [0 0 0] );
            %obj.analysis.fixClosestNode( [l 0], ["uy"], [0] );
            obj.analysis.printProblemInfo();
            obj.plotModel();
        end

    end
end



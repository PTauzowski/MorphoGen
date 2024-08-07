classdef ColumnModel < ModelLinear
    
    methods
        function obj = ColumnModel(sf,l,h,nh,E,nu,P,xp)
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addRectMesh2D( 0, 0, l, h, round(nh*l/h), nh, sf.pattern );
            obj.x=ones(size(obj.mesh.elems,1),1);
            fixedEdgeSelector = Selector( @(x)( abs(x(:,2)) < 0.001 ) );
            loadedEdgeSelector = Selector( @(x)( abs(x(:,2)-h) < 0.001 ) );

            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );            
            %obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true ); 
            obj.analysis = SecondOrderElasticityWeighted( obj.fe, obj.mesh, false );
            %obj.analysis.loadClosestNode(xp, ["ux" "uy"], P );
            obj.analysis.elementLoadLineIntegral( "global", loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + P ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
            %obj.analysis.fixClosestNode( [0 0], ["ux" "uy"], [0 0] );
            %obj.analysis.fixClosestNode( [l 0], ["uy"], [0] );
            obj.analysis.printProblemInfo();
            obj.plotModel();
        end

    end
end



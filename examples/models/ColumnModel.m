classdef ColumnModel < ModelLinear
    
    methods
        function obj = ColumnModel(sf,b,h,l,nh,E,nu,P,xp)
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addRectMesh2D( 0, 0, b, l, round(nh*b/l), nh, sf.pattern );
            obj.x=ones(size(obj.mesh.elems,1),1);
            fixedEdgeSelector = Selector( @(x)( abs(x(:,2)) < 0.001 ) );
            loadedEdgeSelector = Selector( @(x)( abs(x(:,2)-l) < 0.001 ) );

            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.props.h=h;
            obj.fe.setMaterial( material );            
            %obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true ); 
            obj.analysis = SecondOrderElasticityWeighted( obj.fe, obj.mesh, false );
            %obj.analysis.loadClosestNode(xp, ["ux" "uy"], P );
            obj.analysis.elementLoadLineIntegral( "global", loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + P ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
            %obj.analysis.fixNodes( loadedEdgeSelector, ["ux" ] );
            %obj.analysis.fixClosestNode( [0 0], ["ux" "uy"], [0 0] );
            %obj.analysis.fixClosestNode( [b/2 0], ["ux" "uy"], [0 0] );
            obj.analysis.printProblemInfo();
            obj.plotModel();
        end

    end
end



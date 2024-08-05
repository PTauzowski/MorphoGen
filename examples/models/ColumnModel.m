classdef ColumnModel < ModelLinearLoad
    
    methods
        function obj = ColumnModel(sf,l,h,nh,E,nu,P,xp)
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addRectMesh2D( 0, 0, l, h, round(nh*l/h), nh, sf.pattern );
            fixedEdgeSelector = Selector( @(x)( abs(x(:,2)) < 0.001 ) );
            loadedEdgeSelector = Selector( @(x)( abs(x(:,2)-h) < 0.001 ) );

            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );            
            obj.analysis = LinearStability( obj.fe, obj.mesh );           
            %obj.analysis.loadClosestNode(xp, ["ux" "uy"], P );
            obj.analysis.elementLoadLineIntegral( "global", loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + (x(:,1)-l/2).*P ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
            %obj.analysis.fixClosestNode( [0 0], ["ux" "uy"], [0 0] );
            obj.analysis.printProblemInfo();
            obj.plotModel();
        end

        function solve(obj, num_eigenvalues)
            solve(obj, num_eigenvalues);

        end
       
    end
end



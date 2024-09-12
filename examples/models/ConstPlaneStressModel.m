classdef ConstPlaneStressModel< ModelLinear
    properties
        loadedEdgeSelector;
    end
    methods
        function obj = ConstPlaneStressModel(sf,l,res,E,nu,p)

            obj.mesh = Mesh();
            obj.mesh.addRectMesh2D(0, 0, 2*l, l, 2*res, res, sf.pattern);
            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
            fixedEdgeSelector = Selector( @(x)( abs(x(:,1) ) ) < 0.001 );
            obj.loadedEdgeSelector = Selector( @(x)( abs(x(:,1) - 2*l ) < 0.001 ) );
   
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + p ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" ] );
            obj.analysis.fixClosestNode([0 0],["uy"], 0);
            obj.analysis.printProblemInfo();
            obj.x=ones(1,obj.analysis.getTotalElemsNumber());
            obj.result_number=17;
        end

         function setupVariables(obj,E,nu,P)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdgeSelector,  ["ux" "uy"], @(x)( x*0 + P ));
        end
       
    end
end

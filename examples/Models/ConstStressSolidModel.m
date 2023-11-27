classdef ConstStressSolidModel< ModelLinear
    properties
        loadedFaceSelector;
    end
    methods
        function obj = ConstStressSolidModel(sf,l,res,E,nu,pressure)
            hlc=0.4;
            obj.mesh = Mesh();
            obj.mesh.addRectMesh3D( 0, 0, 0, l, l*hlc, l*hlc, res, round(res*hlc), round(res*hlc), sf.localNodes);
            %mesh.addRectMeshTetrahedral3D( '6T', [0 0 0], [2.5*l l l], [2.5*res, res, res] )
            obj.fe = SolidElasticElem( sf, obj.mesh.elems );
       
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material)
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            fixedEdgeSelector = Selector( @(x)( abs(x(:,1) ) < 0.0001 ) );
            obj.loadedFaceSelector = Selector( @(x)( abs(x(:,1) - l) < 0.0001 ) );
            obj.analysis.elementLoadSurfaceIntegral( "global", obj.loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [-pressure 0 0] ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux"] );
            obj.analysis.fixClosestNode([0 0 0],["uy" "uz"], 0);
            
            obj.analysis.printProblemInfo();
            obj.x=ones(1,obj.analysis.getTotalElemsNumber());
            obj.result_number=13;
        end

         function setupVariables(obj,E,nu,pressure)
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.elementLoadSurfaceIntegral( "global", obj.loadedFaceSelector, ["ux" "uy" "uz"], @(x)( x*0 + [-pressure 0 0] ));
        end
       
    end
end

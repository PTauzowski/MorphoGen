classdef CantileverSolidModel < ModelLinear
    properties
        loadedFaceSelector;
    end
    methods
        function obj = CantileverSolidModel(sf,l,res,E,nu,p)
            obj.mesh = Mesh();
            obj.mesh.addRectMesh3D( 0, 0, 0, l, l/2, l/2, 2*res, res, res, sf.localNodes);
            obj.fe = SolidElasticElem( sf, obj.mesh.elems );
    
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material)
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            
            fixedFaceSelector = Selector( @(x)( abs(x(:,1) )<0.001 ) );
            obj.loadedFaceSelector = Selector( @(x)( abs(x(:,1) - l )<0.001 ) );
            obj.analysis.elementLoadSurfaceIntegral( "global", obj.loadedFaceSelector, ["ux" "uy"], @(x)( x(:,1:2)*0 + [0 -p] ));
            obj.analysis.fixNodes( fixedFaceSelector, ["ux" "uy" "uz"] );
            
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


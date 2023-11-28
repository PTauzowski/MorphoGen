classdef LShapeSolidModel < ModelLinear
    
    methods
        function obj = LShapeSolidModel(sf,l,res,E,nu,xp,P)
            obj.mesh = Mesh();
            obj.mesh.addLshape3D( l, 0.4*l, res, sf.localNodes );
            fixedFaceSelector = Selector( @(x)( abs(x(:,3) - l)<0.001 ) );
            loadedFaceSelector = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );
            
            obj.fe = SolidElasticElem( sf, obj.mesh.elems );
            
            obj.fe.plot(obj.mesh.nodes);
            obj.fe.props.h=1;
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material)
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            
            obj.xp=xp;
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], P);
            obj.analysis.fixNodes( fixedFaceSelector, ["ux" "uy" "uz"] );   
            obj.analysis.printProblemInfo();

            obj.x=ones(obj.analysis.getTotalElemsNumber());
            obj.result_number=13;
        end

         function setupVariables(obj,E,nu,P)
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], P);
        end
       
    end
end


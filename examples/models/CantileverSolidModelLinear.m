classdef CantileverSolidModelLinear < ModelLinearLoad
    
    methods
        function obj = CantileverSolidModelLinear(sf,w,res,E,nu,xp,resPt)
            obj.mesh = Mesh();
            obj.mesh.addRectMesh3D( 0,0,0, 8*w, w, 4*w, 8*res,res,4*res, sf.localNodes );
            fixedFaceSelector = Selector( @(x)( abs(x(:,1))<0.001 ) );
            %loadedFaceSelector = Selector( @(x)( abs(x(:,1) - l)<0.001 ) );
            
            obj.fe = SolidElasticElem( sf, obj.mesh.elems );
            
            obj.fe.plot(obj.mesh.nodes);
            obj.fe.props.h=1;
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial(material)
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            
            obj.xp=xp;
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], [1 0 0]);
            obj.analysis.createNextRightHandSideVector();
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], [0 1 0]);
            obj.analysis.createNextRightHandSideVector();
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], [0 0 1]);
            obj.analysis.createNextRightHandSideVector();
            obj.analysis.fixNodes( fixedFaceSelector, ["ux" "uy" "uz"] );   
            obj.analysis.printProblemInfo();
            obj.setResultNode(resPt);
            obj.result_number=13;
            obj.setX(ones(obj.analysis.getTotalElemsNumber(),1));
        end

         function setupVariables(obj,E,nu,P)
            material = SolidMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.loadClosestNode(obj.xp,["ux" "uy" "uz"], P);
         end

         function pstress = computePenalizedHMstress(obj,E,nu,pressure,penalty)
            
            obj.analysis.computeElementResults(obj.x);
            pstress=sum(obj.fe.results.nodal.all(:,14).*obj.fe.results.nodal.all(:,obj.result_number))^(1/penalty);
        end
       
    end
end


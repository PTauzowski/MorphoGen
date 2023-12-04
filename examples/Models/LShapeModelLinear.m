classdef LShapeModelLinear < ModelLinear
    
    methods
        function obj = LShapeModelLinear(sf,l,res,E,nu,xp,P)
            c=0.4;
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addLshape(l,c*l,res,sf.pattern);
            fixedEdgeSelector = Selector( @(x)( abs(x(:,2)-l) ) < 0.001 );

            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );           
            obj.analysis.loadClosestNode(xp,["ux" "uy"], P);
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );    
            obj.analysis.printProblemInfo();
            obj.analysis.initializeResults();
            obj.x=ones(obj.analysis.getTotalElemsNumber(),1);
            obj.result_number=17;
         
        end
       
    end
end


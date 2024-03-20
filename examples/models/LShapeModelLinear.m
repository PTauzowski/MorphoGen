classdef LShapeModelLinear < ModelLinearLoad
    
    methods
        function obj = LShapeModelLinear(sf,l,res,E,nu,xp)
            c=0.4;
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addLshape(l,c*l,res,sf.pattern);
            obj.result_node = obj.mesh.findClosestNode(xp);
            fixedEdgeSelector = Selector( @(x)( abs(x(:,2)-l) ) < 0.001 );

            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, true );
            obj.analysis.loadClosestNode(xp,["ux" "uy"], [1 0]);
            obj.analysis.createNextRightHandSideVector();
            obj.analysis.loadClosestNode(xp,["ux" "uy"], [0 1]);
            obj.analysis.createNextRightHandSideVector();
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );  
            %obj.analysis.fixClosestNode([0,l],["ux" "uy"],[0 0]);
            obj.analysis.printProblemInfo();
            obj.analysis.initializeResults();
            obj.x=ones(obj.analysis.getTotalElemsNumber(),1);
            obj.result_number=17;
            obj.setOneX();
            obj.P0fem=obj.analysis.Pfem;
        end
       
    end
end


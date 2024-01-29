classdef CantileverModelLinear < ModelLinearLoad
    
    methods
        function obj = CantileverModelLinear(sf,l,nh,E,nu,xp)
            c=2;
            obj.xp=xp;
            obj.mesh = Mesh();      
            obj.mesh.addRectMesh2D( 0, 0, l, l/c, round(c*nh), nh, sf.pattern );
            fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) < 0.001 ) );

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
            obj.analysis.printProblemInfo();
            obj.x=ones(1,obj.analysis.getTotalElemsNumber());
            obj.result_number=17;
            obj.result_node=obj.mesh.findClosestNode([l 0]);
            obj.setOneX();
            obj.P0fem=obj.analysis.Pfem;
        end
       
    end
end


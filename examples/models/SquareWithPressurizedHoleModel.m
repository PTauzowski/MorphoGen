classdef SquareWithPressurizedHoleModel < ModelLinear
    properties
        circleSelector;
    end
    methods
        function obj = SquareWithPressurizedHoleModel(sf,x,a,hf,res,E,nu,pressure)
            x0 = x(1);
            y0 = x(2);
            obj.mesh = Mesh();
            obj.mesh.addRectWithHoleMesh2D( a, x0, y0, hf, res, sf.pattern );
            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );
            obj.fe.plot(obj.mesh.nodes);
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh,  false );
            obj.circleSelector = Selector( @(x)( abs(((x(:,1) - x0).^2 + (x(:,2) - y0).^2 ) - (a*hf)^2 )) < 0.01 );
            obj.analysis.elementLoadLineIntegral( "local", obj.circleSelector,  ["ux" "uy"], @(x)( x*0 + [ 0 -pressure ] ));
            obj.analysis.fixClosestNode( [ x0-a y0-a], ["ux" "uy"], [0 0] );
            obj.analysis.fixClosestNode([ x0+a y0-a],"uy", 0 );
            
            obj.analysis.printProblemInfo();
            obj.x=ones(obj.analysis.getTotalElemsNumber());
            obj.result_number=17;
        end

        function setupVariables(obj,E,nu,pressure)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.clearCurrentLoad();
            obj.analysis.elementLoadLineIntegral( "local", obj.circleSelector,  ["ux" "uy"], @(x)( x*0 + [ 0 -pressure ] ));
        end
       
    end
end


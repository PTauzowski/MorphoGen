classdef LameProblemModel < ModelLinear
    properties
        circleSelectorOut;
    end

    methods
        function obj = LameProblemModel(sf,r1,r2,div,E,nu,P)
            x0 = 0;
            y0 = 0;
            obj.mesh = Mesh();
            elems = obj.mesh.addRing2D( x0, y0 , r1, r2, div, round(3*pi*div), sf.pattern );
            %mesh.transformMeshDeg2D( [137 0], -90, [-137 0] );
            obj.fe=PlaneStressElem( sf, elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );

            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
            circleSelectorInt = Selector( @(x)( abs(((x(:,1) - x0).^2 + (x(:,2) - y0).^2 ) - (r1)^2 )) < 0.001 );
            obj.circleSelectorOut = Selector( @(x)( abs(((x(:,1) - x0).^2 + (x(:,2) - y0).^2 ) - (r2)^2 )) < 0.001 );
            obj.analysis.elementLoadLineIntegral( "local", obj.circleSelectorOut,  ["ux" "uy"], @(x)( x*0 + P ));
            obj.analysis.fixClosestNode( [ x0-r1 y0], ["ux" "uy"], [0 0] );
            obj.analysis.fixClosestNode([ x0+r1 y0],"uy", 0 );
            obj.analysis.plotCurrentLoad();
            obj.analysis.plotSupport();           
            obj.analysis.printProblemInfo();
            obj.x=ones(obj.analysis.getTotalElemsNumber(),1);
            obj.result_number=17;
        end

        function setupVariables(obj,E,nu,P)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.elementLoadLineIntegral( "local", obj.circleSelectorOut,  ["ux" "uy"], @(x)( x*0 + P ));
        end
       
    end
end
classdef PressurizedSquareWithHoleModel < ModelLinear
    properties
        loadedEdge1,loadedEdge2,loadedEdge3,loadedEdge4;
    end
    methods
        function obj = PressurizedSquareWithHoleModel(sf,x,a,hf,res,E,nu,pressure)
            x0 = x(1);
            y0 = x(2);
           
            obj.mesh = Mesh();
            elems=obj.mesh.addRectWithHoleMesh2D( 10, x0, y0, hf, res, sf.pattern );
            obj.fe=PlaneStressElem( sf, elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(210000, 0.3);
            obj.fe.setMaterial( material );
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
            
            obj.loadedEdge1 = Selector( @(x)( abs(x(:,1) -(x0 - a) ) < 0.001 ) );
            obj.loadedEdge2 = Selector( @(x)( abs(x(:,1) -(x0 + a) ) < 0.001) );
            obj.loadedEdge3 = Selector( @(x)( abs(x(:,2) -(y0 - a) ) < 0.001 ) );
            obj.loadedEdge4 = Selector( @(x)( abs(x(:,2) -(y0 + a) ) < 0.001 ) );
           % circleSelector = Selector( @(x)( abs(((x(:,1) - x0).^2 + (x(:,2) - y0).^2 ) - (a*hf)^2 )) < 0.01 );
            %problem.plotSelectedNodeNumbers( loadedEdgeSelector );
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge1,  ["ux" "uy"], @(x)( x*0 + [ pressure 0]));
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge2,  ["ux" "uy"], @(x)( x*0 + [-pressure 0]));
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge3, ["ux" "uy"], @(x)(  x*0 + [0 pressure ]));
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge4, ["ux" "uy"], @(x)(  x*0 + [0 -pressure ]));
            %analysis.elementLoadLineIntegral( "local", circleSelector,  ["ux" "uy"], @(x)( x*0 + [ 0 -100 ] ));
            
            obj.analysis.fixClosestNode( [ x0-a y0], ["ux" "uy"], [0 0] );
            obj.analysis.fixClosestNode([ x0+a y0],"uy", 0 );
            
            obj.analysis.printProblemInfo();
            obj.x=ones(obj.analysis.getTotalElemsNumber());
            obj.result_number=17;
        end

        function setupVariables(obj,E,nu,pressure)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );  
            obj.analysis.clearCurrentLoad();
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge1,  ["ux" "uy"], @(x)( x*0 + [ pressure 0]));
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge2,  ["ux" "uy"], @(x)( x*0 + [-pressure 0]));
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge3, ["ux" "uy"], @(x)(  x*0 + [0 pressure ]));
            obj.analysis.elementLoadLineIntegral( "global", obj.loadedEdge4, ["ux" "uy"], @(x)(  x*0 + [0 -pressure ]));
        end
       
    end
end


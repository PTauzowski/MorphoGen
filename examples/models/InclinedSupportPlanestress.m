classdef InclinedSupportPlanestress < ModelLinear
    
    properties
        loadedEdgeSelector
    end

    methods
        function obj = InclinedSupportPlanestress(sf,l,endAngle,res,E,nu,P)
            segmentHeigh=l;
            segmentLength=2*l;
            segmentShortening=segmentHeigh/2*tan(endAngle*pi/180);
            obj.mesh = Mesh();
            obj.mesh.addRectMesh2D(0, -segmentHeigh/2, segmentLength, segmentHeigh, 2*res, res, sf.pattern);
            obj.mesh.transformNodesXY( @(x)([ x(:,1)+2*segmentShortening.*x(:,2)./segmentHeigh.*(x(:,1))./segmentLength  x(:,2) ]) )
            obj.fe=PlaneStressElem( sf, obj.mesh.elems );
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            material.setElasticIzoGrad();
            obj.fe.setMaterial( material );
            
            obj.analysis = LinearElasticityWeighted( obj.fe, obj.mesh, false );
            fixedEdgeSelector = Selector( @(x)( abs(x(:,1)) == 0 ) );
            inclinedSupportEdgeSelector = Selector( @(x)( abs(x(:,1)-(segmentLength+2*segmentShortening.*x(:,2)./segmentHeigh ) ) < 0.001 ) );
            obj.loadedEdgeSelector = Selector( @(x)( abs(x(:,2) - segmentHeigh/2) < 0.001 ) );
            obj.analysis.elementLoadLineIntegral("global",obj.loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + P ));
            obj.analysis.fixNodes( fixedEdgeSelector, ["ux" "uy"] );
            obj.analysis.fixNodes( inclinedSupportEdgeSelector, ["uy"] );
            obj.analysis.setNodalSupportRotations( inclinedSupportEdgeSelector, 90-endAngle);
            obj.analysis.printProblemInfo();
            obj.x=ones(1,obj.analysis.getTotalElemsNumber());
            obj.result_number=17;
        end
        
         function setupVariables(obj,E,nu,P)
            material = PlaneStressMaterial('mat1');
            material.setElasticIzo(E, nu);
            obj.fe.setMaterial( material );   
            obj.analysis.elementLoadLineIntegral("global",obj.loadedEdgeSelector, ["ux" "uy"], @(x)( x*0 + P ));
        end
       
    end
end


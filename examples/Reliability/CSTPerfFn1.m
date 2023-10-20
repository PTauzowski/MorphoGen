classdef CSTPerfFn1 < Function
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        E, F, analysis;
    end
    
    methods
        function obj = CSTPerfFn1(F,E,analysis)
            obj = obj@Function(2, [F*0.0001 E*0.0001])
            obj.F=F;
            obj.E=E;
            obj.FEAnalysis=analysis;
        end
        
        function cfn = computeValue(obj,x)
            material.setElasticIzo(x(1), 0.3);
            material.setElasticIzoGrad();
            obj.analysis.felems{1}.setMaterial( material );
            P=analysis.g
            cfn = obj.fn.computeValue(x) / obj.threshold - 1;
        end

    end
end


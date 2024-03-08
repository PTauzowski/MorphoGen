classdef  loadPerformanceFunctionDisp < Function

    properties
        model;
    end

    methods
        function obj = loadPerformanceFunctionDisp(model)
            obj=obj@Function(2,0.00001)
            obj.model=model;   
        end

        function [g, fi] = computeValue(obj,points)
            g=obj.model.computeLinearDisplacement(points)+0.02;
            fi=[];
        end

    end

end


classdef  loadPerformanceFunctionDisp < Function

    properties
        model;
    end

    methods
        function obj = loadPerformanceFunctionDisp(model,dim)
            obj=obj@Function(dim,0.0001)
            obj.model=model;   
        end

        function [g, fi] = computeValue(obj,points)
            g=obj.model.computeLinearDisplacement(points);
            fi=[];
        end

    end

end


classdef  performanceFunctionPlus < Function

    properties
        g,threshold;
    end

    methods
        
        function obj = performanceFunctionPlus(g)
            obj=obj@Function(g.dim,g.eps);
            obj.g=g;
            obj.resetThreshold();
        end
        function resetThreshold(obj)
            obj.threshold=0;
        end

        function [g, fi] = computeValue(obj,points)
            g=obj.g.computeValue(points)+obj.threshold;
            fi=[];
        end

    end

end


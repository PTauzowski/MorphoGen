classdef  performanceFunctionMinus < Function

    properties
        g,threshold;
    end

    methods
        
        function obj = performanceFunctionMinus(g)
            obj=obj@Function(g.dim,g.eps);
            obj.g=g;
            obj.resetThreshold();
        end
        function resetThreshold(obj)
            obj.threshold=0;
        end

        function [g, fi] = computeValue(obj,points)
            g=obj.threshold-obj.g.computeValue(points);
            fi=[];
        end

    end

end


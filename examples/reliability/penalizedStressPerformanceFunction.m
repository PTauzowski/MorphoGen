classdef  penalizedStressPerformanceFunction < Function

    properties
        model, penalty;
    end

    methods
        function obj = penalizedStressPerformanceFunction(model,dim,penalty)
            obj=obj@Function(dim,0.01)
            obj.model=model;   
            obj.penalty=penalty;
        end

        function [g, fi] = computeValue(obj,points)
            g = obj.model.computePenalizedStress(obj.penalty,points);
            fi=[];
        end

    end

end


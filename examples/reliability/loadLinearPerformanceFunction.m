classdef  loadLinearPerformanceFunction < Function

    properties
        model, g, g0;
    end

    methods
        
        function obj = loadLinearPerformanceFunction(model,g)
            obj=obj@Function(2,0.00001)
            obj.model=model;   
            obj.g=g;
            obj.recomputeBaseFunctions();
        end

        function recomputeBaseFunctions(obj)
            %obj.g0 = obj.g.computeValue(eye(obj.g.dim));
            obj.g0=zeros(obj.g.dim,1);
            points=eye(obj.g.dim);
            for k=1:obj.g.dim
                obj.g0(k) = obj.g.computeValue(points(k,:));
            end

        end

        function g = computeValue(obj,points)
            g=points*obj.g0;         
        end

    end

end


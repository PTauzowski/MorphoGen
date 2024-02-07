classdef  loadPerformanceFunctionDisp < Function

    properties
        model;
    end

    methods
        function obj = loadPerformanceFunctionDisp(model)
            obj=obj@Function(2,0.00001)
            obj.model=model;   
        end

        function g = computeValue(obj,points)
            g=0.013+obj.model.computeLinearDisplacement(points);
        end

    end

end


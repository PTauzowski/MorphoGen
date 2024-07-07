classdef  loadCantileverSolidPerformanceFunctionDisp < Function

    properties
        model;
    end

    methods
        function obj = loadCantileverSolidPerformanceFunctionDisp(model)
            obj=obj@Function(3,0.00001)
            obj.model=model;   
        end

        function g = computeValue(obj,points)
            g=obj.model.computeLinearDisplacement(points);         
        end

    end

end


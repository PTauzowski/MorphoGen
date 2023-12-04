classdef  loadPerformanceFunctionP2 < Function

    properties
        model;
    end

    methods
        function obj = loadPerformanceFunctionP2(model)
            obj=obj@Function(2,0.0001)
            obj.model=model;   
        end

        function g = computeValue(obj,points)
            g=zeros(size(points,1),1);
            for k=1:size(points,1)
                u=obj.model.computeDisplacement(210000,0.3,[points(k,1) points(k,2)]);
                g(k)=u(2)-0.003;
            end            
        end

    end

end


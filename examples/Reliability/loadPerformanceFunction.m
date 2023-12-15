classdef  loadPerformanceFunction < Function

    properties
        model;
    end

    methods
        function obj = loadPerformanceFunction(model)
            obj=obj@Function(2,0.0001)
            obj.model=model;   
        end

        function g = computeValue(obj,points)
            g=zeros(size(points,1),1);
            for k=1:size(points,1)
                u=obj.model.computeDisplacement(210000,0.3,[points(k,1) points(k,2)]);
                g(k)=u(1)+0.017;
            end            
        end

    end

end


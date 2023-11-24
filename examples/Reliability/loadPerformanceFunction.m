classdef  loadPerformanceFunction < Function

    properties
        analysis, material, nelem, nnode, nres,loadedEdgeSelectorX,loadedEdgeSelectorY;
    end

    methods
        function obj = loadPerformanceFunction(model)
            obj.model=model;
           
        end

        function g = computeValue(obj,points)
            g=zeros(size(points,1),1);
            for k=1:size(points,1)
                g(k)=obj.model.computeDisplacement(1,0.3,points(k,:));
            end            
        end

    end

end


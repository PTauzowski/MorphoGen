classdef  loadPenalizedStressPerformanceFunction < Function

    properties
        model;
    end

    methods
        function obj = loadPenalizedStressPerformanceFunction(model)
            obj=obj@Function(2,0.0001)
            obj.model=model;   
        end

        function g = computeValue(obj,points)
            g=zeros(size(points,1),1);
            for k=1:size(points,1)
                ps=obj.model.computePenalizedHMstress(210000,0.3,[points(k,1) points(k,2)],6);
                % g(k)=6.52-ps; % res=20
                g(k)=8.41-ps;
            end            
        end

    end

end


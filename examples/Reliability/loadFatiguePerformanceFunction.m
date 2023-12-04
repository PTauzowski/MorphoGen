classdef  loadFatiguePerformanceFunction < Function

    properties
        model, lcf;
    end

    methods
        function obj = loadFatiguePerformanceFunction(model,fatigueDdata)
            obj=obj@Function(2,0.0001)
            obj.model=model;   
            obj.lcf = LowCycleFatigue(fatigueDdata);
        end

        function g = computeValue(obj,points)
            g=zeros(size(points,1),1);
            for k=1:size(points,1)
                ps=obj.model.computeHMstress(210000,0.3,[points(k,1) points(k,2)]);
                %g(k)=400-sHM;
                g(k)=obj.lcf.nCycles(ps)-10000;
            end            
        end

    end

end


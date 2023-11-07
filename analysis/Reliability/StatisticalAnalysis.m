classdef StatisticalAnalysis < handle
    
    properties
        dim, randVars, perfFn;
    end
    
    methods

        function obj = StatisticalAnalysis( randVars, perfFn )
            obj.randVars = randVars;
            obj.perfFn = perfFn;
            obj.dim = size(obj.randVars,2);
        end
        
        function [x, scatter] = solve(obj,N)
            
            x=zeros(N,obj.dim);
            for k=1:obj.dim
                x(:,k)=obj.randVars{k}.randomize(N);
            end
            scatter=obj.perfFn.computeValue(x);

        end
    end
end


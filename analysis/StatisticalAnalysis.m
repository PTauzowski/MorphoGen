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
        
        function results = solve(obj,N)            
            x=zeros(N,obj.dim);
            for k=1:obj.dim
                x(:,k)=obj.randVars{k}.randomize(N);
            end
            results.scatter=obj.perfFn.computeValue(x);
            results.x=x;
        end
    end
end


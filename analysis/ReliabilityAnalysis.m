classdef ReliabilityAnalysis < handle
    
    properties
        randVars;
        g;
        Pf, beta;
    end
    
    methods(Abstract)
      solve(obj);
    end
    
    methods
        function obj = ReliabilityAnalysis( randVars, g  )
            obj.randVars = randVars;
            obj.g = g;
        end
        
        function dim = getDim(obj)
            dim = size(obj.randVars,2);
        end
        
        function x = generateRandomSapmles(obj,nsamples)
            dim = obj.getDim();
            x=zeros(nsamples,dim);
            for k=1:dim
                x(:,k) = obj.randVars{k}.randomize(nsamples);
            end
        end
        
        function [x, r, fi] = generatePerformanceRandomSapmles(obj,nsamples)
            x=obj.generateRandomSapmles(nsamples);
            [r, fi] = obj.g.computeValue( x );
        end
        
        function obj = plotPerformance(obj,x,y,i)
            figure, hold on;
            plot( x(i), y(i),'.');
            plot( x(not(i)), y(not(i)),'.');
        end
        
        function obj = plotPerformance3D(obj,x,y,z,i)
            figure;
            scatter3(x(i), y(i), z(i), '.' );
            scatter3( x(not(i)), y(not(i)), z(not(i)), '.' );
        end
      
    end                 
end


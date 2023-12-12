classdef ReliabilityAnalysis < handle
    
    properties
        randVars;
        g;
        Pf;
        beta;
        success;
        err_msg;
    end
    
    methods(Abstract)
      solve(obj);
    end
    
    methods
        function obj = ReliabilityAnalysis( randVars, g  )
            obj.randVars = randVars;
            obj.g = g;
            obj.err_msg = "OK";
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
        
        function [x, r, fx] = generatePerformanceRandomSapmles(obj,nsamples)
            dim = obj.getDim();
            x=obj.generateRandomSapmles(nsamples);
            r = zeros(nsamples,1);
            for k=1:nsamples
                value = obj.g.computeValue( x(k,:) );
                if  obj.g.success
                    r(k)=value;
                else
                    fx = [ fx; x(k,:) ];
                end
            end
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


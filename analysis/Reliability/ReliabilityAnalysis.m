classdef ReliabilityAnalysis < handle
    
    properties
        g;
        randVars;
        Pf;
        beta;
    end
    
    methods(Abstract)
      solve(obj);
    end
    
    methods
        function obj = ReliabilityAnalysis( g, randVars )
            obj.g = g;
            obj.randVars = randVars;
        end
        
        function dim = getDim(obj)
            dim = size(obj.randVars,2);
        end
        
        function u = transformToU( obj, x )
            dim = obj.getDim();
            u = zeros( 1, dim );
            for k=1:size(obj.randVars, 2)
               u(k)=obj.randVars(k).toU(x(k)); 
            end
        end
        
        function x = transformFromU( obj, u )
            dim = obj.getDim();
            x = zeros( 1, dim );
            for k=1:size(obj.randVars,2)
               x(k)=obj.randVars(k).fromU(u(k)); 
            end
        end
        
        function x = generateRandomSapmles(obj,nsamples)
            dim = obj.getDim();
            x=zeros(nsamples,dim);
            for k=1:dim
                x(:,k) = obj.randVars{k}.random(nsamples);
            end
        end
        
        function [x, r] = generatePerformanceRandomSapmles(obj,nsamples)
            dim = obj.getDim();
            x=obj.generateRandomSapmles(nsamples);
            r = zeros(nsamples,1);
            for k=1:nsamples
                r(k) = obj.g.computeValue( x(k,:) );
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


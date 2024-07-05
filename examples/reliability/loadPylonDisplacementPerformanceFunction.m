classdef  loadPylonDisplacementPerformanceFunction < Function

    properties
        model, lcf;
    end

    methods
        function obj = loadPylonDisplacementPerformanceFunction(model,fatigueDdata)
            obj=obj@Function(3,0.0001)
            obj.model=model;   
            obj.lcf = LowCycleFatigue(fatigueDdata);
        end

        function g = computeValue(obj,points)
            %g=0.0015-obj.model.computeLinearDisplacement(points); % res=16
            g=obj.model.computeLinearDisplacement(points); %res=12
        end

        function tabNCycles(obj,x0,x1,N)
            dx=(x1-x0)/N;
            plNc=[];
            plX=[];
            for k=1:N
                plNc=[plNc obj.lcf.nCycles(x0+k*dx)];
                plX=[plX x0+k*dx];
            end
            min(plNc)
            max(plNc)
            plot(plX,plNc);
            title('Life cycles');
        end

    end

end

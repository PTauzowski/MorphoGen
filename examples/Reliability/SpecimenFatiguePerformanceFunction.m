classdef  SpecimenFatiguePerformanceFunction < Function

    properties
        model, lcf;
    end

    methods
        function obj = SpecimenFatiguePerformanceFunction(model,fatigueDdata)
            obj=obj@Function(2,0.0001)
            obj.model=model;   
            obj.lcf = LowCycleFatigue(fatigueDdata);
        end

        function g = computeValue(obj,points)
            g=zeros(size(points,1),1);
            for k=1:size(points,1)
                ps=obj.model.computeHMstress(210000,0.3,[points(k,1) points(k,2)]);
                %g(k)=obj.lcf.nCycles(ps)-1000; %2D
                g(k)=obj.lcf.nCycles(ps); %3D
            end            
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


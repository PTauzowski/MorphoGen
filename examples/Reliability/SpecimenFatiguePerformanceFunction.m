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

        function [g, success] = computeValue(obj,points)
            g=zeros(size(points,1),1);
            success=true(size(points,1),1);
            for k=1:size(points,1)
                material = PlaneStressMaterial('mat1');
                material.setElasticIzo(points(k,1), 0.3);
                obj.model.fe.setMaterial( material );   
                obj.lcf.E=points(k,1);
                obj.lcf.sy=points(k,2);
                g(k)=obj.lcf.nCycles(obj.model.computeMaxHMstress());
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


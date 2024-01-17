classdef  chocolatePerformanceFunction < Function
    properties
        height,E,nu,model;
    end
  
    methods
        function obj=chocolatePerformanceFunction(height,E,nu)
            obj=obj@Function(3,0.0001);
            obj.height=height;
            obj.E=E;
            obj.nu=nu;
        end

        function g = computeValue(obj,x)
            g = zeros(size(x,1),1);
            for k=1:size(x,1)
                alGanTh=x(k,1)*obj.height;
                ganTh=(1-x(k,1))*obj.height;
                relRoutndNotchDepth=0.3;
                obj.model = ChocolateModel( ganTh, alGanTh, x(k,3), x(k,2), relRoutndNotchDepth, obj.E, obj.nu);
                obj.model.solveWeighted();  
                %obj.model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
                g(k)=obj.model.computeStressObjective();
               % obj.model.FEAP_Export('chocolateOpti.i',obj.model.mesh,obj.model.ganTh,obj.model.intTh,0.1);
            end
        end

    end

end


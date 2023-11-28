classdef  chocolatePerformanceFunction < Function
    properties
        height,E,nu;
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
                model = ChocolateModel( ganTh, alGanTh, x(k,3), x(k,2), relRoutndNotchDepth, obj.E, obj.nu);
                model.solveWeighted();  
                g(k)=model.computeStressObjective();
            end
        end

    end

end


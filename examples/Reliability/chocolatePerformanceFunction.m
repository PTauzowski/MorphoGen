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

        function createModel(obj,x)
            % x(1) - gan, alGan proportion,  x(2) - relNotchDepth ,    x(3) - notchWidth
            alGanTh=x(1,1)*obj.height;
            ganTh=(1-x(1,1))*obj.height;
            relRoutndNotchDepth=0.3;
            %generateMesh( obj, ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth )
            obj.model = ChocolateModel( ganTh, alGanTh, x(1,3), x(1,2), relRoutndNotchDepth, obj.E, obj.nu);
        end

        function g = computeValue(obj,x)
            g = zeros(size(x,1),1);
            for k=1:size(x,1)
                createModel(obj,x(k,:));
                obj.model.solveWeighted();  
                %obj.model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
                g(k)=obj.model.computeStressObjective();
               % obj.model.FEAP_Export('chocolateOpti.i',obj.model.mesh,obj.model.ganTh,obj.model.intTh,0.1);
            end
        end

         function evaluateValue(obj,x)
                createModel(obj,x(1,:));
                obj.model.plotModel();
                obj.model.solveWeighted();  
                obj.model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
                obj.model.FEAP_Export('chocolateOpti.i',obj.model.mesh,obj.model.ganTh,obj.model.intTh,0.1);
         end

         function fullFactorialBoundsPlot(obj, lb, ub)
                createModel(obj,lb);
                obj.model.plotModel();
                title('Lower bounds');
                createModel(obj,ub);
                figure;
                obj.model.plotModel();
                title('Upper bounds');
         end

    end

end


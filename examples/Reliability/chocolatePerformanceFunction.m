classdef  chocolatePerformanceFunction < Function
    properties
        height,E,nu,alphaT,dT,model,stresses;
    end
  
    methods
        function obj=chocolatePerformanceFunction(height,E,nu,alphaT,dT)
            obj=obj@Function(3,0.0001);
            obj.height=height;
            obj.E=E;
            obj.nu=nu;
            obj.alphaT=alphaT;
            obj.dT=dT;

        end

        function createModel(obj,x)
            % x(1) - gan, alGan proportion,  x(2) - relNotchDepth ,    x(3) - notchWidth
            alGanTh=x(1,1)*obj.height;
            ganTh=(1-x(1,1))*obj.height;
            relRoutndNotchDepth=0.3;
            %generateMesh( obj, ganTh, alGanTh, notchWidth, relNotchDepth, relRoutndNotchDepth )
            obj.model = ChocolateModel( ganTh, alGanTh, x(1,3), x(1,2), relRoutndNotchDepth, obj.E, obj.nu,obj.alphaT,obj.dT);
        end

        function g = computeValue(obj,x)
            g = zeros(size(x,1),1);
            obj.stresses=zeros(size(x,1),4);
            for k=1:size(x,1)
                createModel(obj,x(k,:));
                obj.model.solveWeighted();  
                %obj.model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
                [g(k) obj.stresses(k,1) obj.stresses(k,2) obj.stresses(k,3) obj.stresses(k,4)]=obj.model.computeStressObjective();
               % obj.model.FEAP_Export('chocolateOpti.i',obj.model.mesh,obj.model.ganTh,obj.model.intTh,0.1);
            end
        end

         function evaluateValue(obj,x)
                close all;
                createModel(obj,x(1,:));
                obj.model.plotModel();
                obj.model.solveWeighted();  
                obj.model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.1);
                obj.model.FEAP_Export('chocolateOpti.i',obj.model.mesh,obj.model.ganTh,obj.model.alGanTh,obj.model.intTh,0.1);
                [stressObj, sx1, sy1, sx2, sy2]=obj.model.computeStressObjective();
                fprintf('\nsxx1+syy1=%8.8f, sxx2+syy2=%8.8f, objOpt=%8.8f\n',sx1+sy1,sx2+sy2,stressObj);
         end

         

         function evaluateValue2(obj)
                close all;
                obj.model = ChocolateModel( obj.height*0.67, obj.height*0.33, 4, 0.85, 0.5, obj.E, obj.nu,obj.alphaT,obj.dT);
                obj.model.plotModel();
                obj.model.solveWeighted();  
                obj.model.analysis.plotMaps(["uy" "ux" "sxx" "sxy" "syy" "sHM"],0.05);
                obj.model.FEAP_Export('chocolateOpti.i',obj.model.mesh,obj.model.ganTh,obj.model.alGanTh,obj.model.intTh,0.1);
                [stressObj, sx1, sy1, sx2, sy2]=obj.model.computeStressObjective();
                fprintf('\nsxx1+syy1=%8.8f, sxx2+syy2=%8.8f, objOpt=%8.8f\n',sx1+sy1,sx2+sy2,stressObj);
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


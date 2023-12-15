classdef ReliabilityTaskTuner < handle
    
    properties
        model, topOpt, randVars, g, transform, N, mcDD, mcTop, form, hmv;
        mcDD_res, mcTop_res;
        formDD_res, formTop_res;
        hmvDD_res, hmvTop_res;
    end
    
    methods
        
        function obj = ReliabilityTaskTuner(model, topOpt, randVars, transform, g, N, betat)
            obj.model = model;
            obj.topOpt = topOpt;
            obj.randVars = randVars;
            obj.g=g;
            obj.N=N;
            obj.mcDD  = MonteCarlo(randVars,g,N);
            obj.mcTop = MonteCarlo(randVars,g,N);
            obj.topOpt.is_silent=true;
            obj.hmv = HMV(randVars,g,transform, betat);
            obj.form = FORM(randVars,g,transform);
        end
        
        function obj = tuneMC(obj)
            obj.model.x(:)=1;
            fprintf("\n* Topology optimization with reliability constraints tuner *\n");
            fprintf("Design Domain Monte Carlo sampling\n");
            obj.mcDD_res = obj.mcDD.solve();
            fprintf("Topology optimization iterations\n")
            [~, xopt]  = obj.topOpt.solve();
            obj.model.x=xopt;
            fprintf("Topology Monte Carlo sampling\n");
            obj.mcTop_res = obj.mcTop.solve();
            fprintf("Design domain boundaries: min=%5.7f, max=%5.7f\n",min(obj.mcDD.r),max(obj.mcDD.r));
            fprintf("     Topology boundaries: min=%5.7f, max=%5.7f\n",min(obj.mcTop.r),max(obj.mcTop.r));
        end

         function obj = tuneFORM(obj)
            obj.model.x(:)=1;
            fprintf("\n* FORM with reliability constraints tuner *\n");
            fprintf("Design Domain FORM\n");
            obj.formDD_res = obj.form.solve();
             if size( obj.topOpt.allx, 2 ) > 1
                 obj.model.x=obj.topOpt.allx(:,end);
             else
                 fprintf("Topology optimization iterations\n")
                [~, xopt]  = obj.topOpt.solve();
                obj.model.x=xopt;
             end            
            fprintf("Topology FORMg\n");
            obj.formTop_res = obj.form.solve();
            if obj.formDD_res.success
                fprintf("Design domain FORM: beta=%1.7f\n",obj.formDD_res.beta);
            else
                fprintf("Design domain FORM: not succeed\n");
            end
            if obj.formTop_res.success
                fprintf("     Topology FORM: beta=%1.7f\n",obj.formTop_res.beta);
            else
                fprintf("     Topology FORM: not succeed\n");
            end
         end

         function obj = plotMCs(obj,varNames,objName)
                obj.mcDD.scatterPlots(varNames,['DD' objName])
                obj.mcTop.scatterPlots(varNames,['Top' objName])
         end
    end
end

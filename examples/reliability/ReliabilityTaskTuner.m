classdef ReliabilityTaskTuner < handle
    
    properties
        model, topOpt, randVars, g, transform, N, mcDD, mcTop, mcTopSafe, form, hmv;
        mcDD_res, mcTop_res, mcTopSafe_res;
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
            obj.mcTopSafe = MonteCarlo(randVars,g,N);
            obj.topOpt.is_silent=true;
            obj.hmv = HMV(randVars,g,transform, betat);
            obj.form = FORM(randVars,g,transform);
        end

        function obj = checkTopology(obj,P)
             vfDest=obj.topOpt.Vend;
             obj.topOpt.Vend=0.2;
             obj.model.setupLoad(P);
             [~, ~]  = obj.topOpt.solve();
             obj.topOpt.Vend=vfDest;
        end

        function mpps = checkModality(obj,perf)
            mpps=[];
            P = cellfun(@(x) x.mean, obj.randVars);
            for k=1:obj.g.dim
                Pc=P;
                Pc(k)=Pc(k)+perf*obj.randVars{k}.sd;
                form_res=obj.form.solve(Pc);
                mpps=[mpps; form_res.mpp ];
                Pc=P;
                Pc(k)=Pc(k)-perf*obj.randVars{k}.sd;
                form_res=obj.form.solve(Pc);
                mpps=[mpps; form_res.mpp ];
            end
        end

        function [xmin, xmax] = getBoundsTopologies(obj, vmin, vmax)
             xmin=-1; kmin=-1;
             xmax=-1; kmax=-1;
             plotV=[];
             for k=1:size(obj.topOpt.allx,2)-1
                 vk=sum( obj.topOpt.allx(:,k) )/obj.topOpt.V0;
                 vk1=sum( obj.topOpt.allx(:,k+1) )/obj.topOpt.V0;
                 if vk>=vmin && vk1 < vmin
                     kmin=k;
                     xmin=obj.topOpt.allx(:,k);
                 end
                 if vk>=vmax && vk1 < vmax
                     kmax=k;
                     xmax=obj.topOpt.allx(:,k);
                 end
                 plotV=[plotV vk];
             end
             figure;
             plot(plotV);
             title('Volume fraction evolution');
             if kmin==-1
                 fprintf("\n* ERROR! minimal volume fraction not found *\n");
                 return;
             end
             if kmax==-1
                 fprintf("\n* ERROR! maximal volume fraction not found *\n");
                 return;
             end
             obj.topOpt.setFrame(kmin);
             obj.topOpt.plotCurrentFrame();
             title(["Minimum volume fraction topology" num2str(vmin)]);
             figure;
             obj.topOpt.setFrame(kmax);
             obj.topOpt.plotCurrentFrame();
             title(["Maximum volume fraction topology" num2str(vmax)]);
        end
        
        function obj = fullReliabilityTuning(obj,P,xmin,xmax)
            fprintf("\n* Topology optimization with reliability constraints tuner *\n");
            fprintf("Design Domain Monte Carlo sampling\n");
            obj.g.threshold=0;
            obj.model.setOneX();
            obj.mcDD_res = obj.mcDD.solve();
            obj.mcDD.printResults();
            fprintf("\nMinimum vol stasistics\n");
            obj.model.setX(xmin);
            obj.mcDD_res = obj.mcDD.solve();
            obj.mcDD.printResults();
            fprintf("\nMaximum vol stasistics\n");
            obj.model.setX(xmax);
            obj.mcDD_res = obj.mcDD.solve();
            obj.mcDD.printResults();
            % mpp=obj.hmv.transform.toX(zeros(1,obj.g.dim));
            % obj.model.setupLoad(P);
            % [~, xopt]  = obj.topOpt.solve();
            % obj.model.setX(xopt);
            % fprintf("\nTopology Monte Carlo sampling\n");
        end
        

        function beta = computeBeta(obj,x)
             obj.g.threshold=x;
             form_res = obj.form.solve();
             beta=form_res.beta;
        end


        function threshold = findPerformanceThreshold(obj,m,s,targetBeta)
            backup = obj.g.threshold;
            x1=m-6*s;
            x2=m+6*s;
            b1=obj.computeBeta(x1)-targetBeta;
            b2=obj.computeBeta(x2)-targetBeta;
            while abs(x1-x2)>0.0001*s
                if  b1*b2 < 0 
                    xs=(x1+x2)/2;
                    vs = obj.computeBeta(xs)-targetBeta;
                    if vs*b1 > 0
                        b1=vs;
                        x1=xs;
                    else
                        b2=vs;
                        x2=xs;
                    end
                elseif b1>b2
                    xs=x2+(x2-x1);
                    b2 = obj.computeBeta(xs)-targetBeta;
                    x2 = xs;
                else
                    xs=x1-(x2-x1);
                    b1 = obj.computeBeta(xs)-targetBeta;
                    x1 = xs;
                end
            end
            obj.g.threshold=xs;
            form_res = obj.form.solve();
            form_res.beta
            form_res.mpp
            threshold=xs;
            % obj.g.threshold=0;
            % fn_g = @(x)( (obj.computeBeta(x)-targetBeta)^2 );
            % threshold = fmincon(fn_g,m,[],[],[],[],xmin,xmax);
            % obj.g.threshold = backup;

        end

        function obj = tuneMC(obj)
            fprintf("\n* Topology optimization with reliability constraints tuner *\n");
            fprintf("Design Domain Monte Carlo sampling\n");
            obj.mcDD_res = obj.mcDD.solve();
            obj.mcDD.printStats();
            fprintf("\nTopology optimization iterations\n");
            mpp=obj.hmv.transform.toX(zeros(1,obj.g.dim));
            obj.model.setupLoad(mpp);
            %obj.g.setupLoad(mpp);
            [~, xopt]  = obj.topOpt.solve();
            obj.model.setX(xopt);
            fprintf("\nTopology Monte Carlo sampling\n");
            obj.mcTop_res = obj.mcTop.solve();
            obj.mcTop.printStats();
        end

         function obj = tuneFORM(obj)
            obj.model.setOneX()
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
            obj.hmvTop_res = obj.hmv.solve();
            obj.model.x=xopt;
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
                obj.mcDD.scatterPlots(varNames,['DD' objName]);
                obj.mcTop.scatterPlots(varNames,['Top' objName]);
         end
        
         function obj = plotMC3D(obj,v1,v2,varNames,objName)
                obj.mcDD.scatterPlot3D(v1,v2,varNames,['DD' objName]);
                obj.mcTop.scatterPlot3D(v1,v2,varNames,['Top' objName]);
         end

         function tabPf(obj,x0,min,max,n)
             x=min;
             dx=(max-min)/n;
             for k=1:n
                 obj.g.threshold=x+k*dx;
                 %mc_res = obj.mcDD.solve();
                 form_res=obj.form.solve(x0);
                 %fprintf("x=%1.7f, mean=%1.7f, sd=%1.7f, beta=%1.7f\n",x+k*dx, mc_res.mv, mc_res.sd, mc_res.beta);
                 fprintf("x=%1.7f, beta=%1.7f, mpp:",x+k*dx, form_res.beta);
                 for l=1:size(form_res.mpp,2)
                    fprintf("x(%1d)=%3.4f ",l,form_res.mpp(l));    
                 end
               fprintf("\n");
             end
         end

         
    end
end


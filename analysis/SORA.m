classdef SORA < handle
     properties
        model, topOpt, randVars, g, transform, form, hmv, mpps, betat, baseName, allXi, x0;
    end
    
    methods
        function obj = SORA(baseName, model, topOpt, randVars, g, transform, betat )
            obj.baseName=baseName;
            obj.model=model;
            obj.topOpt=topOpt;
            obj.randVars=randVars;
            obj.form = FORM(randVars,g,transform);
            obj.hmv = HMV(randVars,g,transform,betat);  
            obj.g=g;
            obj.transform=transform;
            obj.topOpt.is_silent=true;
            obj.betat=betat;
            obj.x0=zeros(1,g.dim);
            for k=1:g.dim
                obj.x0(k)=randVars{k}.mean;
            end
        end

        function obj = tune(obj,N)
            obj.model.setOneX();
            obj.g.threshold=0;

            fprintf("\nDesign domain reliability \n");
            mc=MonteCarlo(obj.randVars,obj.g,N);
            mc_resDD=mc.solve();
            form_resDD=obj.form.solve(obj.x0); 
            hmv_resDD=obj.hmv.solve(obj.x0);
            mc.printResults("Design domain",mc_resDD);
            obj.hmv.printResults("Design domain",hmv_resDD);
            obj.form.printResults("Design domain",form_resDD);
            mc.printStats();

            fprintf("\nMean value topology\n");
            obj.model.setOneX();
            obj.model.setupLoad(obj.x0);
            obj.topOpt.solve();

            fprintf("Mean Value reliability\n");
            mc_resMV=mc.solve();
            form_resMV=obj.form.solve(obj.x0);
            hmv_resMV=obj.hmv.solve(obj.x0);
            obj.g.threshold=-(mc_resMV.mv-mc_resMV.sd);
            mc.printResults("Mean value",mc_resMV);
            obj.hmv.printResults("Mean value",hmv_resMV);
            obj.form.printResults("Mean value",form_resMV);
            title('Mean value  topology');

            figure;
            fprintf("\nFirst mpp topology\n");
            obj.model.setOneX();
            obj.model.setupLoad(hmv_resMV.mpp);
            obj.topOpt.solve();

            fprintf("First mpp reliability\n");
            mc_resFirst=mc.solve();
            form_resFirst=obj.form.solve(obj.x0);
            hmv_resFirst=obj.hmv.solve(obj.x0);
            mc.printResults("First mpp",mc_resFirst);
            obj.hmv.printResults("First mpp",hmv_resFirst);
            obj.form.printResults("First mpp",form_resFirst);
            title('First mpp topology');
            obj.tabReliability();

            figure;
            fprintf("\nSecond mpp topology\n");
            obj.model.setOneX();
            obj.model.setupLoad(hmv_resFirst.mpp);
            obj.topOpt.Vend=0.1;
            obj.topOpt.solve();

            fprintf("Second mpp reliability\n");
            mc_resSec=mc.solve();
            form_resSec=obj.form.solve(obj.x0);
            hmv_resSec=obj.hmv.solve(obj.x0);
            mc.printResults("Second  mpp",mc_resSec);
            obj.hmv.printResults("Second  mpp",hmv_resSec);
            obj.form.printResults("Second  mpp",form_resSec);
            title('Second mpp topology');
            fprintf("\n");
        end
    
        function plotG(xp,dx)
            [X,Y]=meshgrid(1:0.5:10,1:20);
            Z = sin(X) + cos(Y);
            surf(X,Y,Z)

            scatter3(obj.x(obj.r>0,1),obj.x(obj.r>0,2),obj.r(obj.r>0),'MarkerEdgeColor',[0 .6 .0],'Marker','.');
            scatter3(obj.x(obj.r<=0,1),obj.x(obj.r<=0,2),obj.r(obj.r<=0),'filled','MarkerEdgeColor',[0.5 0 .5],'Marker','o');
        end

        function tabReliability(obj)
             plBetaPred=[];
             plBetaFORM=[];
             plVolFr=[];
             plG=[];
             plG0=[];
             for k=1:size(obj.topOpt.allx,2)
               %[g0, gmpp, mpp, mpp0_pred, beta_pred] = obj.computePerformance(k);
                form_resMV=obj.form.solve(obj.x0);
                hmv_resMV=obj.hmv.solve(obj.x0);
                g0=hmv_resMV.g0;
                gmpp=hmv_resMV.g;
                mpp=hmv_resMV.mpp;
                mpp0_pred=hmv_resMV.mpp0_pred;
                beta_pred=hmv_resMV.beta_pred;
               obj.model.setX(obj.topOpt.allx(:,k));
               form_res=obj.form.solve(obj.x0);
               obj.topOpt.setFrame(k);
               volFr=obj.topOpt.computeVolumeFraction();
               fprintf("Iter :%3d  VolFr=%3.1f, G0=%5.7f, G=%5.7f, G0-G=%5.7f, BetaPred=%1.4f, BetaFORM=%1.4f\n",k,volFr,g0,gmpp,g0-gmpp,beta_pred,form_res.beta);
               plBetaPred=[plBetaPred beta_pred];
               plBetaFORM=[plBetaFORM form_res.beta];
               plVolFr=[plVolFr volFr];
               plG=[plG gmpp];
               plG0=[plG0 g0];
             end
                figure;
             plot(plVolFr);
             title('Volume fraction evolution ');
             figure;
             plot(plBetaPred);
             title('HMV predicted beta evolution ');
             figure;
             plot(plBetaFORM);
             title('FORM beta evolution ');
             figure, hold on;
             plot(plG);
             plot(plG0);
             title('Performance functions evolution ');
             figure;
             plot(smooth(plG));
             title('Smoothed performance function evolution ');
             figure;
             plot(plVolFr,plBetaFORM);
             title('FORM beta vs volFr ');
        end
       
    end
end


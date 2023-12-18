classdef SORA < handle
    
    properties
        model, topOpt, randVars, form, hmv, mpps, betat;
    end
    
    methods
        function obj = SORA(model, topOpt, randVars, g, transform, betat )
            obj.model=model;
            obj.topOpt=topOpt;
            obj.randVars=randVars;
            obj.form = FORM(randVars,g,transform);
            obj.hmv = HMV(randVars,g,transform,betat);  
            obj.topOpt.is_silent=true;
            obj.betat=betat;
            figure;
        end

     
        
        function results = solveX(obj)
            mpp=obj.hmv.transform.toX([0 0]);
            obj.model.setupLoad(mpp);
            inProgress=true;
            iter=1;
            while inProgress
                if iter==2
                    figure, hold on;
                end
                obj.topOpt.allx=[];
                [~, xopt] = obj.topOpt.solve();  
                
                fr_res=obj.findBetaFrame();
                obj.topOpt.setFrame( fr_res.frame );    
                obj.topOpt.plotCurrentFrame();
                conv=norm(fr_res.mpp-mpp);
                %obj.model.x=fr_res.x;
                %form_res=obj.form.solve();
                fprintf("\nconv=%1.5f, frame=%3d/%3d, beta_form=%5.7f, g=%5.7f, beta_pred=%5.7f, vol=%5.7f, ",conv,fr_res.frame,fr_res.lastframe,fr_res.beta_pred,fr_res.g,fr_res.beta_pred,obj.topOpt.computeVolumeFraction());
                fprintf("mpp: ");
                for k=1:size(mpp,2)
                   fprintf("x(%1d)=%3.4f ",k,mpp(k));
                end   
                if conv<0.001
                    fr = obj.form.solve();
                    results.mpp=mpp;
                    results.xopt=fr_res.x;
                    results.beta_form=fr.beta;
                    fr_res=obj.findBetaFrame();
                    return;
                end
                mpp=fr_res.mpp;
                obj.model.setupLoad(fr_res.mpp);
                %obj.setMeans(mpp);
                iter=iter+1;
            end
        end

        function [g0, gmpp, mpp, mpp0_pred, beta_pred] = computePerformance(obj, iter)
             x=obj.model.x;
             if iter== size(obj.topOpt.allx,2)
                obj.model.x=obj.topOpt.allx(:,iter);
             else
                obj.model.x=obj.topOpt.allx(:,iter);
                res_hmv = obj.hmv.solve();
                g1=res_hmv.g;
                obj.model.x=obj.topOpt.allx(:,iter+1);
                res_hmv = obj.hmv.solve();
                g2=res_hmv.g;
                alpha=abs(g1)/(abs(g2)+g1);
                obj.model.x=alpha*obj.topOpt.allx(:,iter)+(1-alpha)*obj.topOpt.allx(:,iter+1);
             end
             res_hmv = obj.hmv.solve();
             g0=res_hmv.g0;
             gmpp=res_hmv.g;
             beta_pred=res_hmv.beta_pred;
             mpp=res_hmv.mpp;
             %umpp=obj.form.transform.toU(mpp);
             %n=umpp/norm(umpp);
             mpp0_pred=res_hmv.mpp0_pred;
             obj.model.x=x;
        end

        function tabMultiMpp(obj)
            form_res=obj.form.solve();
            cf=0.25;
            plBetaFORM=[];
            plVolFr=[];
            for k=1:40
              umpp=obj.form.transform.toU(form_res.mpp);
              n=umpp/norm(umpp);
              mpp0_pred=obj.form.transform.toX(cf*k*obj.betat*n);
              obj.model.setupLoad(mpp0_pred);
              obj.topOpt.allx=[];
              [~, xopt] = obj.topOpt.solve();  
              fr_res=obj.findBetaFrame();
              obj.topOpt.setFrame( fr_res.frame );    
              obj.topOpt.plotCurrentFrame();
              fprintf("\niter :%3d, cf=%2.5f, frame=%3d/%3d, beta_form=%5.7f, g=%5.7f, beta_pred=%5.7f, vol=%5.7f",k,cf*k,fr_res.frame,fr_res.lastframe,fr_res.beta_pred,fr_res.g,fr_res.beta_pred,obj.topOpt.computeVolumeFraction());
                fprintf(" mpp: ");
                for k=1:size(mpp0_pred,2)
                   fprintf("x(%1d)=%3.4f ",k,mpp0_pred(k));
                end 
                volFr=obj.topOpt.computeVolumeFraction();
                plBetaFORM=[plBetaFORM fr_res.beta_pred];
                plVolFr=[plVolFr volFr];
            end
            figure;
            plot(plVolFr,plBetaFORM);
            title('FORM beta vs volFr');
        end

        function tabReliability(obj)
             plBetaPred=[];
             plBetaFORM=[];
             plVolFr=[];
             plG=[];
             for k=1:size(obj.topOpt.allx,2)
               [g0, gmpp, mpp, mpp0_pred, beta_pred] = obj.computePerformance(k);
               obj.model.x=obj.topOpt.allx(:,k);
               form_res=obj.form.solve();
               obj.topOpt.setFrame(k);
               volFr=obj.topOpt.computeVolumeFraction();
               fprintf("Iter :%3d  VolFr=%3.1f, G=%5.7f, BetaPred=%1.4f, BetaFORM=%1.4f\n",k,volFr,gmpp,beta_pred,form_res.beta);
               plBetaPred=[plBetaPred beta_pred];
               plBetaFORM=[plBetaFORM form_res.beta];
               plVolFr=[plVolFr volFr];
               plG=[plG gmpp];
             end
                figure;
             plot(plVolFr);
             title('Volume fraction evolution');
             figure;
             plot(plBetaPred);
             title('HMV predicted beta evolution');
             figure;
             plot(plBetaFORM);
             title('FORM beta evolution');
             figure;
             plot(plG);
             title('Performance function evolution');
             figure;
             plot(smooth(plG));
             title('Smoothed performance function evolution');
             figure;
             plot(plVolFr,plBetaFORM);
             title('FORM beta vs volFr');
        end

        function limitReliability(obj)
                % obj.topOpt.allx=[];
                % mpp=obj.hmv.transform.toX([0 0]);
                % obj.model.setupLoad(mpp);
                % obj.topOpt.solve();
               [g0, gmpp, mpp, mpp0_pred, beta_pred] = obj.computePerformance(1);
               volFr=sum(obj.topOpt.allx(:,1))/size(obj.topOpt.allx,1)*100;
               fprintf("\nIter :%3d  VolFr=%3.1f, G0=%5.7f,  Gmpp=%5.7f, Beta=%1.4f, ",1,volFr,g0,gmpp,beta_pred);
               for k=1:size(mpp,2)
                  fprintf("x(%1d)=%3.4f ",k,mpp(k));
               end
               [g0, gmpp, mpp, mpp0_pred, beta_pred]  = obj.computePerformance(size(obj.topOpt.allx,2));
               volFr=sum(obj.topOpt.allx(:,size(obj.topOpt.allx,2)))/size(obj.topOpt.allx,1)*100;
               fprintf("\nIter :%3d  VolFr=%3.1f, G0=%5.7f,  Gmpp=%5.7f, Beta=%1.4f, ",1,volFr,g0,gmpp,beta_pred);
               for k=1:size(mpp,2)
                  fprintf("x(%1d)=%3.4f ",k,mpp(k));
               end
               fprintf("\n");
        end

        function results = findBetaFrame(obj)
            results.lastframe=size(obj.topOpt.allx,2);
            ia=1; 
            [ga, gmpp, mpp, mpp0_pred, beta_pred] = obj.computePerformance(ia);
            if gmpp < 0
               results.frame=1;
               results.x=obj.topOpt.allx(:,1);
               results.mpp=mpp;
               results.g=gmpp;
               results.beta_pred=beta_pred;
               results.mpp0_pred=mpp0_pred;
               return;
            end
            ib=size(obj.topOpt.allx,2);
            [gb, gmpp, mpp, mpp0_pred, beta_pred] =obj.computePerformance(ib);
            if gmpp >= 0
               results.frame=ib;
               results.x=obj.topOpt.allx(:,ib);
               results.mpp=mpp;
               results.g=gmpp;
               results.beta_pred=beta_pred;
               results.mpp0_pred=mpp0_pred;
               return;
            end
            while ib-ia>1
                i=round((ib+ia)/2);
                [gb, v, mpp, mpp0_pred, beta_pred]=obj.computePerformance(i);
                if v > 0
                    ia=i;
                    va=v;
                else
                    ib=i;
                    vb=v;
                end
                results.mpp=mpp;
            end
            [gb, v, mpp, mpp0_pred, beta_pred]=obj.computePerformance(ia);
            [gb1, v1, mpp1, mpp0_pred1, beta_pred1]=obj.computePerformance(ia+1);
            results.frame=ia;
            results.x=obj.topOpt.allx(:,ia);
            results.g=v;
            results.beta_pred=beta_pred;
            results.mpp0_pred=mpp0_pred;
        end

        function ci = findZeroGold(obj)
            ia=1; 
            va=obj.computePerformance(ia);
            if va < 0
               ci=1;
               return;
            end
            ib=size(obj.topOpt.allx,2);
            vb=obj.computePerformance(ib);
            if vb > 0
               ci=-1;
               return;
            end
            while ib-ia>1
                i=round((ib+ia)/2);
                v=obj.computePerformance(i);
                if v>0
                    ia=i;
                    va=v;
                else
                    ib=i;
                    vb=v;
                end
            end
            ci=ia;
        end

        function obj=setMeans(obj,m)
            for k=1:size(obj.randVars,2)
                obj.randVars{k}.setMean(m(k));
            end
        end
    end
end


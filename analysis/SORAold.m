classdef SORAold < handle
    
    properties
        model, topOpt, randVars, g, transform, form, hmv, mc, mc_res, mpps, betat, baseName, allXi, x0, xtop;
    end
    
    methods
        function obj = SORAold(baseName, model, topOpt, randVars, g, transform, betat )
            obj.baseName=baseName;
            obj.model=model;
            obj.topOpt=topOpt;
            obj.randVars=randVars;
            obj.form = FORM(randVars,g,transform);
            obj.hmv = HMV(randVars,g,transform,betat);  
            obj.mc = MonteCarlo(randVars, g, 1000000000);
            obj.g=g;
            obj.transform=transform;
            obj.topOpt.is_silent=true;
            obj.betat=betat;
            obj.x0=zeros(1,g.dim);
            for k=1:g.dim
                obj.x0(k)=randVars{k}.mean;
            end
        end

     
        
        function results = solveX(obj,destVol)
            mpp=obj.hmv.transform.toX(zeros(1,obj.g.dim));
            mpp2=mpp;
            obj.model.setupLoad(mpp);
            inProgress=true;
            iter=1;
            convp=10000;
            while inProgress
                if iter==2
                    figure, hold on;
                end
                obj.topOpt.allx=[];
                [~, xopt] = obj.topOpt.solve();  
                fr_res=obj.findBetaFrame();
                %fr_res=obj.checkLastFrame();
                obj.topOpt.setFrame( fr_res.frame );    
                obj.topOpt.plotCurrentFrame();
                if iter==1
                     obj.model.setX(obj.topOpt.allx(:,end));
                    form_res=obj.form.solve(obj.x0);
                    if form_res.success
                        title([obj.baseName 'initial topology, volfr=' num2str(sum(fr_res.x)/size(fr_res.x,1)) 'beta FORM = ' num2str(form_res.beta)]);
                        fprintf('\nDeterministic topology beta FORM=%1.5f',form_res.beta );
                    else
                        title([obj.baseName 'initial topology, volfr=' num2str(sum(fr_res.x)/size(fr_res.x,1))]);
                    end
                    obj.allXi=obj.topOpt.allx;
                    savefig([obj.baseName '_topt.fig']);

                end
                conv=norm(obj.transform.toU(fr_res.mpp)-obj.transform.toU(mpp));
                conv2=norm(obj.transform.toU(fr_res.mpp)-obj.transform.toU(mpp2));
                fprintf("\nconv=%1.5f, conv2=%1.5f, frame=%3d/%3d, beta_form=%5.7f, g=%5.7f, beta_pred=%5.7f, vol=%5.7f, ",conv,conv2,fr_res.frame,fr_res.lastframe,fr_res.beta_pred,fr_res.g,fr_res.beta_pred,obj.topOpt.computeVolumeFraction());
                fprintf("mpp: ");
                for k=1:size(mpp,2)
                   fprintf("x(%1d)=%3.4f ",k,mpp(k));
                end   
                if conv<0.01 || abs(conv-convp)<0.0001 || conv2 < 0.01 || iter==5
                    destFrame=obj.findVolFrame(destVol);
                    obj.topOpt.setFrame(destFrame);
                    obj.xtop=obj.topOpt.allx(:,destFrame);
                    form_res=obj.form.solve(obj.x0);
                    obj.mc_res=obj.mc.solve();
                    obj.topOpt.plotCurrentFrame();
                    title([obj.baseName 'Safe topology, volfr = ' num2str(sum(fr_res.x)/size(fr_res.x,1))]);
                    savefig([obj.baseName '_safe.fig']);
                    %obj.model.setX(obj.topOpt.allx(:,end));
                    obj.model.setX(obj.topOpt.allx(:,destFrame));
                    if form_res.success
                        fprintf('\nProbablistic topology beta FORM=%1.5f',form_res.beta );
                    else
                        fprintf('\nProbablistic topology beta FORM=not succeed');
                    end
                    
                    fprintf("\n\n");
                    obj.allXi=obj.topOpt.allx;
                    return;
                end
                mpp2=mpp;
                mpp=fr_res.mpp;
                obj.model.setupLoad(fr_res.mpp);
                iter=iter+1;
                convp=conv;
            end
            
        end

        function results = solveMean(obj)
            mpp0=obj.hmv.transform.toX(zeros(1,obj.g.dim));
            obj.model.setupLoad(mpp);
            inProgress=true;
            iter=1;
            while inProgress
                obj.topOpt.allx=[];
                [~, xopt] = obj.topOpt.solve();  
            end
        end

        function [g0, gmpp, mpp, mpp0_pred, beta_pred] = computePerformance(obj, iter)
             x=obj.topOpt.x;
             if iter== size(obj.topOpt.allx,2)
                obj.model.setX(obj.topOpt.allx(:,iter));
             else
                obj.model.setX(obj.topOpt.allx(:,iter));
                % res_hmv = obj.hmv.solve(obj.x0);
                % g1=res_hmv.g;
                % obj.model.x=obj.topOpt.allx(:,iter+1);
                % res_hmv = obj.hmv.solve(obj.x0);
                % g2=res_hmv.g;
                % alpha=abs(g1)/(abs(g2)+g1);
                % obj.model.x=alpha*obj.topOpt.allx(:,iter)+(1-alpha)*obj.topOpt.allx(:,iter+1);
             end
             res_hmv = obj.hmv.solve(obj.x0);
             if res_hmv.success
                g0=res_hmv.g0;
                gmpp=res_hmv.g;
                beta_pred=res_hmv.beta_pred;
                mpp=res_hmv.mpp;
                %umpp=obj.form.transform.toU(mpp);
                %n=umpp/norm(umpp);
                mpp0_pred=res_hmv.mpp0_pred;
                obj.model.setX(x);
             else
                 a=1;
             end
        end

        function tabMultiMpp(obj)
            form_res=obj.form.solve(obj.x0);
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
            title('FORM beta vs volFr ');
        end

        function tabReliability(obj)
             plBetaPred=[];
             plBetaFORM=[];
             plVolFr=[];
             plG=[];
             plG0=[];
             for k=1:size(obj.topOpt.allx,2)
               [g0, gmpp, mpp, mpp0_pred, beta_pred] = obj.computePerformance(k);
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

        function frame = findVolFrame(obj,vf)
            frame=0;
             for k=1:size(obj.topOpt.allx,2)
                 frame=k;
                 if vf>sum(obj.topOpt.allx(:,k))/obj.topOpt.V0
                     return;
                 end
             end
        end

        function limitReliability(obj)
               obj.topOpt.allx=[];
               mpp=obj.hmv.transform.toX(zeros(1,obj.g.dim));               
               obj.model.setupLoad(mpp);
               obj.topOpt.solve();
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
            ib=size(obj.topOpt.allx,2);
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
            [g01, g1, mpp1, mpp0_pred1, beta_pred1]=obj.computePerformance(ia);
            [g02, g2, mpp2, mpp0_pred2, beta_pred2]=obj.computePerformance(ia+1);
            x2=g1/(g1-g2);
            x1=-g2/(g1-g2);
            results.frame=ia;
            results.x=obj.topOpt.allx(:,ia)*x1+obj.topOpt.allx(:,ia+1)*x2;
            results.g=x1*g1+x2*g2;
            results.beta_pred=beta_pred1*x1+beta_pred2*x2;
            results.mpp0_pred=mpp0_pred1*x1+mpp0_pred2*x2;
        end

       
        function results = checkLastFrame(obj)
            results.lastframe=size(obj.topOpt.allx,2);
            i=size(obj.topOpt.allx,2);
            [g0, v, mpp, mpp0_pred, beta_pred]=obj.computePerformance(i);
            results.frame=i;
            results.x=obj.topOpt.allx(:,end);
            results.g=v;
            results.beta_pred=beta_pred;
            results.mpp0_pred=mpp0_pred;
            results.mpp=mpp;
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

        function checkTuning(obj)
            plBetaTop=[];
            plBetaSafe=[];
            plVolTop=[];
            plVolSafe=[];

            zeroLoad=obj.transform.toX( zeros(1,obj.g.dim) );
%            obj.model.setOneX();
            fprintf("\n* FORM with reliability constraints tuner *\n");
            fprintf("Design Domain FORM\n");
            formDD_res = obj.form.solve(obj.x0);
            obj.model.setupLoad(zeroLoad);
            if formDD_res.success
                fprintf("Design domain  FORM: beta=%1.7f\n",formDD_res.beta);
            else
                fprintf("Design domain FORM: not succeed\n");
             end
             fprintf("Topology optimization ...\n");
             obj.allXi=obj.topOpt.allx;
             obj.topOpt.allx=[];
             obj.topOpt.solve();           
             fr_res=obj.findBetaFrame();
             obj.topOpt.setFrame( fr_res.frame );    
             obj.topOpt.plotCurrentFrame();
             obj.model.setX(fr_res.x);
             obj.topOpt.x=obj.model.getX();
             formTop_res = obj.form.solve(obj.x0);
             fprintf("Topology optimization FORM\n");
             if formTop_res.success
                fprintf("     Topology, frame=%3d/%3d, FORM: beta=%1.7f, g=%5.7f, beta_pred=%5.7f, vol=%5.7f\n",fr_res.frame,fr_res.lastframe,formTop_res.beta, fr_res.g,fr_res.beta_pred,obj.topOpt.computeVolumeFraction());
             else
                fprintf("     Topology FORM: not succeed\n");
             end

             fprintf("Topology optimization ...\n");
             obj.model.setupLoad(fr_res.mpp);
             obj.topOpt.allx=[];
             obj.topOpt.solve();
            
             fr_res=obj.findBetaFrame();
             obj.topOpt.setFrame( fr_res.frame );    
             obj.topOpt.plotCurrentFrame();
             obj.model.setX(fr_res.x);
             obj.topOpt.x=obj.model.getX();
             formTop_res = obj.form.solve(obj.x0);
             fprintf("Safe topology optimization FORM\n");
             if formTop_res.success
                 fprintf("Safe Topology, frame=%3d/%3d, FORM: beta=%1.7f, g=%5.7f, beta_pred=%5.7f, vol=%5.7f\n",fr_res.frame,fr_res.lastframe,formTop_res.beta, fr_res.g,fr_res.beta_pred,obj.topOpt.computeVolumeFraction());
             else
                fprintf("Safe topology FORM: not succeed\n");
             end
             figure, hold on;
             plot(plBetaTop);
             plot(plBetaSafe);
             figure, hold on;
             plot(plVolTop);
             plot(plVolSafe);
        end
    end
end


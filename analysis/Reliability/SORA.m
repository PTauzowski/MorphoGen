classdef SORA < handle
    
    properties
        model, topOpt, form, hmv, mpps, betat;
    end
    
    methods
        function obj = SORA(model, topOpt, randomVariables, g, transform, betat )
            obj.model=model;
            obj.topOpt=topOpt;
            obj.form = FORM(randomVariables,g,transform);
            obj.hmv = HMV(randomVariables,g,transform,betat);  
            obj.topOpt.is_silent=true;
            obj.betat=betat;
            % if isfile('_allXopt.mat')
            %      load('_allXopt.mat','allx');
            %      obj.topOpt.allx=allx;
            %      obj.topOpt.x=allx(:,end);
            %      obj.topOpt.setFrame( size(allx,2) );
            %      obj.topOpt.plotCurrentFrame();
            % else
                  [objF, xopt] = obj.topOpt.solve();
                  allx=obj.topOpt.allx;
                  %save('_allXopt.mat','allx');
            %end
            model.x=obj.topOpt.allx(:,size(allx,2));
            figure;
        end

        function results = solve(obj)
            inProgress=true;
            hmv_res = obj.hmv.solve();
            mpp=0.*hmv_res.mpp;
            mpp=obj.hmv.transform.toX(mpp);
            while inProgress
                obj.model.setupLoad(mpp);
                obj.topOpt.allx=[];
                [objF, xopt] = obj.topOpt.solve(); 
                hmv_res = obj.hmv.solve();
                mpp=hmv_res.mpp;
                obj.mpps=[obj.mpps; mpp];
                fprintf("mppX: ");
                for k=1:size(mpp,2)
                   fprintf("x(%1d)=%3.4f ",k,mpp(k));
                end                
                ci = obj.findBetaFrame();
                [g0, gmpp, mpp1, beta_pred] = obj.computePerformance(ci);
                obj.topOpt.setFrame( ci );                
                obj.topOpt.plotCurrentFrame();
                cn=norm(mpp1-mpp);
                fprintf("conv=%1.5f, g=%5.7f, vol=%5.7f\n",cn,hmv_res.g,obj.topOpt.computeVolumeFraction());
                if cn < 0.001
                    results.mpp=mpp;
                    obj.model.setupLoad(mpp);
                    [g0, gmpp, mpp, beta_pred] = obj.computePerformance(ci);
                    results.mpp=mpp;
                    results.g0=g0;
                    results.beta_pred=beta_pred;
                    results.gmpp=gmpp;
                    results.form_res=obj.form.solve();
                    return;
                end         
                mpp=mpp;
            end
        end
        
        function [betaDet, xRel, volDet, volRel] = solveOld(obj)
            inProgress=true;
            %[ ~, mpp, betaDet ] = obj.reliability.solve();
            [gz, mpp, betaDet] = obj.computePerformance(size(obj.topOpt.allx,2));
            volDet=obj.topOpt.computeVolumeFraction();
            fprintf("* SORA *\n",gz,betaDet,volDet);
            fprintf("g=%5.7f, betaDet=%5.7f, volDet=%5.7f",gz,betaDet,volDet);
            figure;
            mpps=[];
            while inProgress
                mppX=obj.hmv.transform.toX(mpp);
                obj.mpps=[obj.mpps; mppX];
                obj.model.setupLoad(mppX);
                fprintf("\nmppX: ");
                for k=1:size(mppX,2)
                    fprintf("x(%1d)=%3.4f ",k,mppX(k));
                end
                obj.topOpt.allx=[];
                [objF, xopt] = obj.topOpt.solve();         
                zeroiter = obj.findZeroG();
                if zeroiter==-1
                    fprintf("The initial structure is safe. No need to optimize\n");
                    return;
                end
                [gz, mpp1, beta] = obj.computePerformance(zeroiter);
                obj.topOpt.setFrame( zeroiter );                
                obj.topOpt.plotCurrentFrame();
                
                conv=norm(mpp-mpp1);
                fprintf("conv=%3.4f, beta=%1.5f",conv,beta);
                if conv<0.001 
                    inProgress=false;
                    volRel=obj.topOpt.computeVolumeFraction();
                    xRel=obj.topOpt.allx(:,zeroiter);
                else
                    mpp=mpp1;
                end
               
            end
        end

        function [g0, gmpp, mpp, beta_pred] = computePerformance(obj, iter)
             obj.model.x=obj.topOpt.allx(:,iter);
             res_hmv = obj.hmv.solve();
             g0=obj.hmv.g.computeValue( obj.hmv.transform.toX(res_hmv.mpp*0) );
             gmpp=res_hmv.g;
             beta_pred=res_hmv.beta_pred;
             mpp=res_hmv.mpp;
        end

        function tabReliability(obj)
             for k=1:size(obj.topOpt.allx,2)
               [g0, gmpp, mpp, beta_pred] = obj.computePerformance(k);
               volFr=sum(obj.topOpt.allx(:,k))/size(obj.topOpt.allx,1)*100;
               fprintf("Iter :%3d  VolFr=%3.1f, G=%5.7f, Beta=%1.4f\n",k,volFr,gmpp,beta_pred);
             end
        end

        function limitReliability(obj)
               [g0, gmpp, mpp, beta_pred] = obj.computePerformance(1);
               volFr=sum(obj.topOpt.allx(:,1))/size(obj.topOpt.allx,1)*100;
               fprintf("\nIter :%3d  VolFr=%3.1f, G0=%5.7f,  Gmpp=%5.7f, Beta=%1.4f, ",1,volFr,g0,gmpp,beta_pred);
               for k=1:size(mpp,2)
                  fprintf("x(%1d)=%3.4f ",k,mpp(k));
               end
               [g0, gmpp, mpp, beta_pred]  = obj.computePerformance(size(obj.topOpt.allx,2));
               volFr=sum(obj.topOpt.allx(:,size(obj.topOpt.allx,2)))/size(obj.topOpt.allx,1)*100;
               fprintf("\nIter :%3d  VolFr=%3.1f, G0=%5.7f,  Gmpp=%5.7f, Beta=%1.4f, ",1,volFr,g0,gmpp,beta_pred);
               for k=1:size(mpp,2)
                  fprintf("x(%1d)=%3.4f ",k,mpp(k));
               end
               fprintf("\n");
        end

        function ci = findBetaFrame(obj)
            ia=1; 
            [ga, gmpp, mpp, beta_pred] = obj.computePerformance(ia);
            if ga < 0
               ci=1;
               return;
            end
            ib=size(obj.topOpt.allx,2);
            [gb, gmpp, mpp, beta_pred] =obj.computePerformance(ib);
            if gb > 0
               ci=ib;
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
           

    end
end


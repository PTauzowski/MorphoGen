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
            if isfile('_allXopt.mat')
                 load('_allXopt.mat','allx');
                 obj.topOpt.allx=allx;
                 obj.topOpt.x=allx(:,end);
                 obj.topOpt.setFrame( size(allx,2) );
                 obj.topOpt.plotCurrentFrame();
            else
                  [objF, xopt] = obj.topOpt.solve();
                  allx=obj.topOpt.allx;
                  save('_allXopt.mat','allx');
            end
            model.x=obj.topOpt.allx(:,size(allx,2));
        end

        function [betaDet, xRel, volDet, volRel] = solve(obj)
            inProgress=true;
            [ ~, mpp_hmv, betar ] = obj.hmv.solve();
            mpp=0.*mpp_hmv;
            while inProgress
                [ ~, mpp_hmv, betar ] = obj.hmv.solve();
                dmpp=mpp_hmv-mpp;
                mpp1=mpp+dmpp*(betar-obj.betat);
                mppX=obj.hmv.transform.toX(mpp1);
                obj.mpps=[obj.mpps; mppX];
                obj.model.setupLoad(mppX);
                fprintf("\nmppX: ");
                for k=1:size(mppX,2)
                   fprintf("x(%1d)=%3.4f ",k,mppX(k));
                end
                obj.topOpt.allx=[];
                [objF, xopt] = obj.topOpt.solve(); 
                ci = obj.findBetaFrame();
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

        function [g, mpp, beta] = computePerformance(obj, iter)
             obj.model.x=obj.topOpt.allx(:,iter);
             [ ~, mpp, beta ] = obj.hmv.solve();
             g=obj.hmv.g.computeValue( obj.hmv.transform.toX(mpp) );
        end

        function tabReliability(obj)
             for k=1:size(obj.topOpt.allx,2)
               [g, ~, betar] = obj.computePerformance(k);
               volFr=sum(obj.topOpt.allx(:,k))/size(obj.topOpt.allx,1)*100;
               fprintf("Iter :%3d  VolFr=%3.1f, G=%5.7f, Beta=%1.4f\n",k,volFr,g,betar);
             end
        end

        function limitReliability(obj)
               [g, mpp, betar] = obj.computePerformance(1);
               volFr=sum(obj.topOpt.allx(:,1))/size(obj.topOpt.allx,1)*100;
               fprintf("\nIter :%3d  VolFr=%3.1f, G=%5.7f, Beta=%1.4f, ",1,volFr,g,betar);
               mppX=obj.hmv.transform.toX(mpp);
               for k=1:size(mppX,2)
                  fprintf("x(%1d)=%3.4f ",k,mppX(k));
               end
               [g, mpp, betar] = obj.computePerformance(size(obj.topOpt.allx,2));
               volFr=sum(obj.topOpt.allx(:,size(obj.topOpt.allx,2)))/size(obj.topOpt.allx,1)*100;
               fprintf("\nIter :%3d  VolFr=%3.1f, G=%5.7f, Beta=%1.4f, ",size(obj.topOpt.allx,2),volFr,g,betar);
               mppX=obj.hmv.transform.toX(mpp);
               for k=1:size(mppX,2)
                  fprintf("x(%1d)=%3.4f ",k,mppX(k));
               end
               fprintf("\n");
        end

        function ci = findBetaFrame(obj)
            ia=1; 
            [ga, mppa, betaa] = obj.computePerformance(ia);
            if ga < 0
               ci=1;
               return;
            end
            ib=size(obj.topOpt.allx,2);
            [gb, mppb, betab] =obj.computePerformance(ib);
            if gb > 0
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


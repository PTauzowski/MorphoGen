classdef SORA
    
    properties
        model, topOpt, reliability;
    end
    
    methods
        function obj = SORA(model, topOpt, reliability)
            obj.model=model;
            obj.topOpt=topOpt;
            obj.reliability=reliability;
            obj.topOpt.is_silent=true;
        end
        
        function [xDet, betaDet, xRel, volDet, volRel] = solve(obj)
            inProgress=true;
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
            %[ ~, mpp, betaDet ] = obj.reliability.solve();
            [gz, mpp, betaDet] = obj.computePerformance(size(allx,2));
            betaDet
            volDet=obj.topOpt.computeVolumeFraction()
            xDet=obj.topOpt.allx(:,end);
            figure;
            while inProgress
                mppX=obj.reliability.transform.toX(mpp);
                obj.model.setupLoad(mppX);
                [objF, xopt] = obj.topOpt.solve();
                % for k=1:size(obj.topOpt.allx,2)
                %    [g, ~, betar] = obj.computePerformance(k);
                %    fprintf("Iter :%3d  G=%5.7f, Beta=%1.4f\n",k,g,betar);
                % end
                % obj.reliability.transform.toX(mpp)
                zeroiter = obj.findZeroG();
                [gz, mpp1, beta] = obj.computePerformance(zeroiter);
                obj.topOpt.setFrame( zeroiter );
                obj.topOpt.plotCurrentFrame();
                % gz1 = obj.computePerformance(zeroiter+1);
                % fprintf("ZeroIterter :%3d,   Gz=%5.7f,  Gz1=%5.7f\n",zeroiter,gz,gz1);
                conv=norm(mpp-mpp1);
                if conv<0.001 
                    inProgress=false;
                    volRel=obj.topOpt.computeVolumeFraction();
                    xRel=obj.topOpt.allx(:,zeroiter);
                else
                    mpp=mpp1;
                    fprintf("conv=%5.7f, beta=%1.5f\n",conv,beta);
                end
            
            end
        end

        function [g, mpp, beta] = computePerformance(obj, iter)
             obj.model.x=obj.topOpt.allx(:,iter);
             [ ~, mpp, beta ] = obj.reliability.solve();
             g=obj.reliability.g.computeValue( obj.reliability.transform.toX(mpp) );
        end

        function ci = findZeroG(obj)
            ia=1; 
            va=obj.computePerformance(ia);
            ib=size(obj.topOpt.allx,2);
            vb=obj.computePerformance(ib);
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


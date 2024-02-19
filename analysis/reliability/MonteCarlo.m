classdef MonteCarlo < ReliabilityAnalysis
        
    properties
        nsamples;
        x;
        r;
        fi;
    end
    
    methods
        function obj = MonteCarlo(randVars, g, nsamples)
            obj = obj@ReliabilityAnalysis(randVars,g);
            obj.nsamples=nsamples;
        end
        
        function results = solve(obj)
            obj.x=obj.generateRandomSapmles(obj.nsamples);
            [obj.x, obj.r, obj.fi] = obj.generatePerformanceRandomSapmles(obj.nsamples);
            results.Pf = max(size(find(obj.r<0)))/(obj.nsamples-size(obj.fi,2));
            results.beta=-norminv(results.Pf);
            obj.Pf=results.Pf;
            obj.beta=results.beta;
        end

        function scatterPlots(obj,varNames,objName)
            for k=1:size(obj.x,2)
                figure, hold on;
                scatter(obj.x(obj.r>0,k),obj.r(obj.r>0),'MarkerEdgeColor',[0 .6 .0],'Marker','.');
                scatter(obj.x(obj.r<=0,k),obj.r(obj.r<=0),'MarkerEdgeColor',[0.5 0 0.5],'Marker','o');
                title(strjoin([varNames(k) " x "  objName]));
            end

            if size(obj.x,2)==2
                figure, hold on;
                scatter3(obj.x(obj.r>0,1),obj.x(obj.r>0,2),obj.r(obj.r>0),'MarkerEdgeColor',[0 .6 .0],'Marker','.');
                scatter3(obj.x(obj.r<=0,1),obj.x(obj.r<=0,2),obj.r(obj.r<=0),'filled','MarkerEdgeColor',[0.5 0 .5],'Marker','o');
                title(strjoin([varNames(1) " x " varNames(2) " x " objName]));
            end
        end
        function printResults(obj)
            mv=mean(obj.r);
            sd=std(obj.r);
            fprintf(" mean=%5.7f, sigma=%5.7f\n",mv,sd);
            fprintf(" min=%5.7f, max=%5.7f\n",min(obj.r),max(obj.r));
            fprintf(" mean -/+  sigma=%5.7f, %5.7f\n",mv-sd,mv+sd);
            fprintf(" mean -/+ 2sigma=%5.7f, %5.7f\n",mv-2*sd,mv+2*sd);
            fprintf(" mean -/+ 3sigma=%5.7f, %5.7f\n",mv-3*sd,mv+3*sd);
            fprintf(" mean -/+ 6sigma=%5.7f, %5.7f\n",mv-6*sd,mv+6*sd);
            fprintf(" Pf=%5.7f, beta=%5.7f\n",obj.Pf, obj.beta);
        end
        
    end
end


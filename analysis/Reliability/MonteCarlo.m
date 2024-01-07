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
    end
end


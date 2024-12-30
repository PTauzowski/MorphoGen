classdef StressIntensityMultiAnalysisTopologyOptimization < StressIntensityTopologyOptimizationVol
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FEAnalyses, plLambda, plVol, lastStableFrame;
    end
    
    methods
        function obj = StressIntensityMultiAnalysisTopologyOptimization(Rmin,FEAnalyses,maxais,penal,volFr,is_const)
            obj=obj@StressIntensityTopologyOptimizationVol(Rmin,FEAnalyses(1),maxais,penal,volFr,is_const);
            obj.FEAnalyses=FEAnalyses;
            obj.maxais=maxais;
            obj.penal=penal;
            obj.x(:)=1;
            obj.xfat=obj.x;
            obj.pnormfat=100;
            obj.plLambda = [];
            obj.plVol = [];
        end
        
        function printIterationInfo(obj)
            fprintf('%5i ',obj.iteration);
            fprintf('Vrel=%2.1f ',round(sum( obj.x )/obj.V0*1000)/10);
            fprintf('lambda=%5.3g ', obj.FEAnalysis.lambda);
            fprintf('\n');
            obj.plLambda = [ obj.plLambda abs(obj.FEAnalysis.lambda) ];
            obj.plVol = [ obj.plVol round(sum( obj.x )/obj.V0*1000)/10 ];
            if abs(obj.FEAnalysis.lambda)>=1
                obj.lastStableFrame=obj.iteration;
            end
        end

        function ais = computeAverageIntensities(obj)
            nFEAnalyses=size(obj.FEAnalyses,2);
            sais=zeros(obj.FEAnalysis.getTotalElemsNumber(),nFEAnalyses);
            ais=zeros(obj.FEAnalysis.getTotalElemsNumber(),1);
            for k=1:nFEAnalyses
                obj.qnodal = obj.FEAnalyses(k).solveWeighted((obj.x).^obj.penal);
                obj.FEAnalyses(k).computeElementResults(obj.x.^obj.penal);
                for i=1:size(obj.FEAnalysis.felems,2)
                   hmIndex=find(obj.FEAnalyses(k).felems{i}.results.names == "sHM");
                   for j=1:size(obj.FEAnalyses(k).felems{i}.elems,1)
                        %ais(obj.elem_inds{i}(j),k) = mean( obj.linearElasticProblem.felems{i}.results.GPvalues(hmIndex,j,:) );
                        sais(obj.elem_inds{i}(j),k) = mean( obj.FEAnalyses(k).felems{i}.results.nodal.all(obj.FEAnalyses(k).felems{i}.elems(j,:),hmIndex) );
                   end
                end
                obj.maxstress = [ obj.maxstress max(sais(:,k)) ];
                ais = max(ais,  sais(:,k) / max(sais(:,k)) );
            end
        end

    end
end


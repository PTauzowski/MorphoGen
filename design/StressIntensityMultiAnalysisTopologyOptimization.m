classdef StressIntensityMultiAnalysisTopologyOptimization < StressIntensityTopologyOptimization
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FEAnalyses;
    end
    
    methods
        function obj = StressIntensityMultiAnalysisTopologyOptimization(numberOfConstraints,Rmin,FEAnalyses,maxais,penal,is_const)
            obj=obj@StressIntensityTopologyOptimization(numberOfConstraints,Rmin,FEAnalyses{1},is_const);
            obj.FEAnalyses=FEAnalyses;
            obj.maxais=maxais;
            obj.penal=penal;
            obj.x(:)=1;
            obj.xfat=obj.x;
            obj.pnormfat=100;
        end

        function ais = computeAverageIntensities(obj)
            nFEAnalyses=size(obj.FEAnalyses,2);
            sais=zeros(obj.FEAnalysis.getTotalElemsNumber(),nFEAnalyses);
            ais=zeros(obj.FEAnalysis.getTotalElemsNumber(),1);
            for k=1:nFEAnalyses
                obj.qnodal = obj.nFEAnalyses{k}.solveWeighted((obj.x).^obj.penal);
                obj.FEAnalyses{k}.computeElementResults(obj.x.^obj.penal);
                for i=1:size(obj.FEAnalysis.felems,2)
                   hmIndex=find(obj.FEAnalyses{k}.felems{i}.results.names == "sHM");
                   for j=1:size(obj.FEAnalyses{k}.felems{i}.elems,1)
                        %ais(obj.elem_inds{i}(j),k) = mean( obj.linearElasticProblem.felems{i}.results.GPvalues(hmIndex,j,:) );
                        sais(obj.elem_inds{i}(j),k) = mean( obj.FEAnalyses{k}.felems{i}.results.nodal.all(obj.FEAnalyses{k}.felems{i}.elems(j,:),hmIndex) );
                   end
                end
                obj.maxstress = [ obj.maxstress max(sais(:,k)) ];
                ais = max(ais,  sais(:,k) / max(ais) );
            end
        end

    end
end


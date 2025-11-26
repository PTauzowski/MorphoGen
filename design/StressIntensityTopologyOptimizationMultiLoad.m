classdef StressIntensityTopologyOptimizationMultiLoad < StressIntensityTopologyOptimizationVol
    
    properties
        reduction_fn;
    end

    methods
            
        function obj = StressIntensityTopologyOptimizationMultiLoad(Rmin,FEAnalysis,maxais,penal,reduction_fn,vol,is_const)
            obj=obj@StressIntensityTopologyOptimizationVol(Rmin,FEAnalysis,maxais,penal,vol,is_const);
            obj.reduction_fn=reduction_fn;
        end
                       
        

    end
end


classdef StressIntensityMultiMaxTopologyOptimization < StressIntensityMultiTopologyOptimization
    methods
        function obj = StressIntensityMultiMaxTopologyOptimization(Rmin, FEAnalyses, maxais, penal, volFr, is_const)
            obj = obj@StressIntensityMultiTopologyOptimization( ...
                Rmin, FEAnalyses, maxais, penal, volFr, is_const, "max");
        end
    end
end

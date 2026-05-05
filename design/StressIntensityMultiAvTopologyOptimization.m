classdef StressIntensityMultiAvTopologyOptimization < StressIntensityMultiTopologyOptimization
    methods
        function obj = StressIntensityMultiAvTopologyOptimization(Rmin, FEAnalyses, maxais, penal, volFr, is_const)
            obj = obj@StressIntensityMultiTopologyOptimization( ...
                Rmin, FEAnalyses, maxais, penal, volFr, is_const, "sum");
        end
    end
end

classdef SIMP_MMA_TopologyOptimizationElasticCompliance < SIMP_MMA_TopologyOptimizationElasticComplianceBase
    methods
        function obj = SIMP_MMA_TopologyOptimizationElasticCompliance(Rmin, problem, penal, VolConstr, is_const)
            obj = obj@SIMP_MMA_TopologyOptimizationElasticComplianceBase( ...
                Rmin, problem, penal, VolConstr, is_const, "linear");
        end
    end
end

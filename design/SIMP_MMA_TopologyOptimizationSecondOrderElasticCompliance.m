classdef SIMP_MMA_TopologyOptimizationSecondOrderElasticCompliance < SIMP_MMA_TopologyOptimizationElasticComplianceBase
    methods
        function obj = SIMP_MMA_TopologyOptimizationSecondOrderElasticCompliance(Rmin, problem, penal, VolConstr, is_const)
            obj = obj@SIMP_MMA_TopologyOptimizationElasticComplianceBase( ...
                Rmin, problem, penal, VolConstr, is_const, "secondOrder");
        end
    end
end

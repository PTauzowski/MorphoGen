classdef ReliabilityConstrainedTopologyOptimization < StressIntensityTopologyOptimization
    
    properties
         reliability, pfEnd;
    end
    
    methods
        function obj = ReliabilityConstrainedTopologyOptimization(Rmin,linearElasticProblem,maxais,penal,reliability,pfEnd,is_const)
            obj=obj@StressIntensityTopologyOptimization(1,Rmin,linearElasticProblem,maxais,penal,is_const)
            obj.pfEnd=pfEnd;
            obj.reliability=reliability;
        end
                           
        function of = computeObjectiveFunction(obj)
            of = sum( obj.x );
            obj.FobjValue=of;
        end

        function dc = computeInequalityConstraints(obj,x)
            [ Pf, ~, ~ ] = obj.reliability.solve();
            dc = Pf - obj.PfEnd;
        end

        function dc = computeEqualityConstraints(obj)
        end
        
        function printIterationInfo(obj)
            fprintf('%5i ',obj.iteration);
            fprintf('Vrel=%2.1f ',round(sum( obj.x )/obj.V0*1000)/10);
            fprintf('\n');
        end
        
    end
end


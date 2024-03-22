classdef StressIntensityTopologyOptimizationG < StressIntensityTopologyOptimization
    
    properties
         g;
    end
    
    methods
        function obj = StressIntensityTopologyOptimizationG(Rmin,linearElasticProblem,maxais,penal,g,is_const)
            obj=obj@StressIntensityTopologyOptimization(1,Rmin,linearElasticProblem,maxais,penal,is_const)
            obj.g=g;
        end
                           
        function of = computeObjectiveFunction(obj)
            of = sum( obj.x );
            obj.FobjValue=of;
        end

        function dc = computeInequalityConstraints(obj,x)
            dc = obj.g.computeValue( x );
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


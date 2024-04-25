classdef FatigueConstrainedTopologyOptimization < StressIntensityTopologyOptimization
    
    properties
         fatigue, Nt;
    end
    
    methods
        function obj = FatigueConstrainedTopologyOptimization(Rmin,linearElasticProblem,maxais,penal,fatigue,Nt,is_const)
            obj=obj@StressIntensityTopologyOptimization(1,Rmin,linearElasticProblem,maxais,penal,is_const)
            obj.Nt=Nt;
            obj.fatigue=fatigue;
        end
                           
        function of = computeObjectiveFunction(obj)
            of = sum( obj.x );
            obj.FobjValue=of;
        end

        function dc = computeInequalityConstraints(obj,x)
            Nc = obj.fatigue.nCycles(obj.maxstress(end));
            dc = Nc - obj.Nt;
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


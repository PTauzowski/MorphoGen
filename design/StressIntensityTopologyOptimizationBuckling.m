classdef StressIntensityTopologyOptimizationBuckling < StressIntensityTopologyOptimization
    
    properties
         Vend,plLambda,plVol,lastStableFrame;
    end
    
    methods
        function obj = StressIntensityTopologyOptimizationBuckling(Rmin,linearElasticProblem,maxais,penal,Vend,is_const)
            obj=obj@StressIntensityTopologyOptimization(1,Rmin,linearElasticProblem,maxais,penal,is_const)
            obj.Vend=Vend;
            obj.plLambda=[];
            obj.plVol=[];
        end
                           
        function of = computeObjectiveFunction(obj)
            of = sum( obj.x );
            obj.FobjValue=of;
        end

        function dc = computeInequalityConstraints(obj,x)
            dc = sum( x )/obj.V0 - obj.Vend;
        end

        function dc = computeEqualityConstraints(obj)
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
        
    end
end


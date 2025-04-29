classdef StressIntensityTopologyOptimizationDynamicBuckling < StressIntensityTopologyOptimization
    
    properties
         Vend,plLambda,plVol,plOmegas,lastStableFrame;
    end
    
    methods
        function obj = StressIntensityTopologyOptimizationDynamicBuckling(Rmin,linearElasticProblem,maxais,penal,Vend,is_const)
            obj=obj@StressIntensityTopologyOptimization(1,Rmin,linearElasticProblem,maxais,penal,is_const)
            obj.Vend=Vend;
            obj.plLambda=[];
            obj.plVol=[];
            obj.plOmegas=[];
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
            fprintf('omega(1)=%5.3g ', obj.FEAnalysis.omegas(1));
            fprintf('omega(2)=%5.3g ', obj.FEAnalysis.omegas(2));
            fprintf('omega(3)=%5.3g ', obj.FEAnalysis.omegas(3));
            fprintf('omega(4)=%5.3g ', obj.FEAnalysis.omegas(4));
            fprintf('\n');
            if abs(obj.FEAnalysis.omegas(1))>10000
                obj.plLambda = [ obj.plLambda abs( obj.FEAnalysis.lambda) ];
                obj.plOmegas = [ obj.plOmegas abs( obj.FEAnalysis.omegas) ];
                obj.plVol = [ obj.plVol round(sum( obj.x )/obj.V0*1000)/10 ];
            end
            
            if abs(obj.FEAnalysis.lambda)>=1
                obj.lastStableFrame=obj.iteration;
            end
        end
        
    end
end


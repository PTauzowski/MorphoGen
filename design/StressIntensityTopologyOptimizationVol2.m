classdef StressIntensityTopologyOptimizationVol2 < StressIntensityTopologyOptimization
    
    properties
         Vend, plotFobj, plotCn, cnodes
    end
    
    methods
        function obj = StressIntensityTopologyOptimizationVol2(Rmin,linearElasticProblem,maxais,penal,Vend,is_const)
            obj=obj@StressIntensityTopologyOptimization(1,Rmin,linearElasticProblem,maxais,penal,is_const)
            obj.Vend=Vend;
            obj.plotFobj=[];
            obj.plotCn=[];
        end
                           
        function of = computeObjectiveFunction(obj)
            of = sum( obj.x );
            obj.FobjValue=of;

        end

        function dc = computeInequalityConstraints(obj,x)
            dc = sum( x )/obj.V0 - obj.Vend;
            obj.plotCn = [ obj.plotCn max(obj.FEAnalysis.qnodal(obj.cnodes,2) ) ];
            obj.plotFobj = [ obj.plotFobj sum( obj.x/obj.V0*100.00 ) ];
            obj.allx=[obj.allx x ];
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


classdef SIMP_MMA_TopologyOptimization < GradientBasedTopologyOptimization  
    
    properties
        penal; 
        a0,ai,ci,di,low,upp,xold1,xold2,change;
    end

   
    methods
        
        function obj = SIMP_MMA_TopologyOptimization(numberOfConstraints,Rmin,FEproblem,penal,is_const)
            obj=obj@GradientBasedTopologyOptimization(numberOfConstraints,Rmin,FEproblem,is_const);
            obj.penal=penal;
            obj.xold1(1:obj.totalFENumber,1) = 0; 
            obj.xold2(1:obj.totalFENumber,1) = 0; 
            obj.a0=1;
            obj.ai(1:numberOfConstraints,1)=0;
            obj.ci(1:numberOfConstraints,1)=1000;
            obj.di(1:numberOfConstraints,1)=0;
            obj.low=obj.ai;
            obj.upp=obj.ai;
            obj.change=1;
        end

        function computeObjectiveFunction(obj, x)

        end

        function rc = isNotFinished(obj)
            rc = obj.change >= 0.001;
        end

        function updateDesign(obj)
            obj.computeObjectiveFunctonWithGradient(obj.x);
            obj.computeConstraintsAndGradient(obj.x);
            obj.gradFobjValue = obj.filteringByMAtrix( obj.gradFobjValue );
            [xmma,~,~,~,~,~,~,~,~,obj.low,obj.upp] = mmasub2( ...
                    size(obj.constrValues,1), ...
                    obj.totalFENumber, ...
                    obj.iteration, ...
                    obj.x, ...
                    obj.xmin, obj.xmax, ...
                    obj.xold1, obj.xold2, ...
                    obj.FobjValue, obj.gradFobjValue, 0*obj.gradFobjValue, ...
                    obj.constrValues,obj.gradConstrValues,0*obj.gradConstrValues, ...
                    obj.low,obj.upp,obj.a0,obj.ai,obj.ci,obj.di);
                if obj.iteration > 1
                    obj.xold2 = obj.xold1;
                end
                obj.xold1=obj.x;
                obj.x=xmma;
                obj.change = max(max(abs(obj.x-obj.xold1)));      
        end
        
    end
end


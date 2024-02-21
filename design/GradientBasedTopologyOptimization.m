classdef GradientBasedTopologyOptimization < TopologyOptimization
    
    properties
        gradFobjValue, gradConstrValues;
    end

    methods

        function obj = GradientBasedTopologyOptimization(numberOfConstraints,Rmin,FEproblem,is_const)
            obj = obj@TopologyOptimization(numberOfConstraints,Rmin,FEproblem,is_const)
            obj.gradFobjValue = zeros(size(obj.x,1),1);
            obj.gradConstrValues=zeros(numberOfConstraints,size(obj.x,1));
        end

        function computeObjectiveFunctionAndGradient(obj,x)
            obj.computeObjectiveFunctionGradient(x);
            obj.computeConstraintsGradient(x);
        end

        function computeConstraintsAndGradient(obj,x)
            obj.computeConstraints(x);
            obj.computeConstraintsGradient(x);
        end

    end
end


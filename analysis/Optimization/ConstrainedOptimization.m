classdef (Abstract) ConstrainedOptimization < Optimization

    properties
        numberOfConstraints, constrValues;
    end

    methods
        function obj = ConstrainedOptimization(numberOfDesignVariables,numberOfConstraints)
            obj=obj@Optimization(numberOfDesignVariables)
            obj.numberOfConstraints;
            obj.constrValues=zeros(numberOfConstraints,1);
        end       
    end

end


classdef (Abstract) Optimization < handle
    
    properties
        numberOfDesignVariables, FobjValue, x, xmin, xmax, iteration;
    end

    methods
        function obj = Optimization(numberOfDesignVariables)
            obj.numberOfDesignVariables = numberOfDesignVariables;
            obj.x=zeros(numberOfDesignVariables,1);
            obj.xmin=zeros(numberOfDesignVariables,1);
            obj.xmax=zeros(numberOfDesignVariables,1);
            obj.iteration=1;
        end

    end

    methods (Abstract)
        computeObjectiveFunction(obj,x);
        isNotFinished(obj);
        solve(obj);
    end

end


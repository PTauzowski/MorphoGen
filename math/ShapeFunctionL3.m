classdef ShapeFunctionL3 < ShapeFunctions
    methods
        function obj = ShapeFunctionL3() 
            obj=obj@ShapeFunctions(1)
            obj.localNodes = [-1; 0; 1];
            obj.vertices=[1 2 3]';
            obj.pattern = [0 1 2]';
        end
        function value = computeValue( ~, xi )
            value = [ xi.*(xi-1)/2 (1-xi.*xi) xi.*(xi+1)/2 ];
        end
        function grad = computeGradient( ~, xi )
            grad = [ (2*xi-1)/2  -2*xi (1+2*xi)/2 ];
        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationRectangular(2,1);
        end
    end
end 

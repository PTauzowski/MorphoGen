classdef ShapeFunctionL2 < ShapeFunctions
    methods
        function obj = ShapeFunctionL2() 
            obj=obj@ShapeFunctions(1);
            obj.localNodes = [-1; 1];
            obj.vertices=[1 2]';
            obj.pattern = [0 1]';
        end
        function value = computeValue( ~, xi )
            np = size( xi, 1 );
            value = [ (1-xi)/2 (1+xi)/2 ];
        end
        function grad = computeGradient( ~, xi )
            grad = repmat([ -1/2 1/2 ],size(xi,1),1);
        end
        
        function integrator = createIntegrator(varargin)
            integrator = IntegrationRectangular(1,1);
        end
        
    end
end 

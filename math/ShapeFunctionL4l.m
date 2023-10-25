classdef ShapeFunctionL4l < ShapeFunctions
    methods
        function obj = ShapeFunctionL4l() 
            obj=obj@ShapeFunctions(1)
            obj.localNodes = [-1; -1/3; 1/3; 1];
            obj.vertices=[1 4]';
            obj.pattern = [0 1 2 3]';
        end
        function value = computeValue( ~, x )
            value = [ -(x-1.0).*(3.0*x-1.0).*(3.0*x+1.0)/16.0...
                       9.0*(x-1.0).*(x+1.0).*(3.0*x-1.0)/16.0 ...
                      -9.0*(x-1.0).*(x+1.0).*(3.0*x+1.0)/16.0 ...
                      (x+1.0).*(3.0*x-1.0).*(3.0*x+1.0)/16.0 ];
         end
        function grad = computeGradient( ~, xi )
            x = xi;
            xx = xi.*xi;
            grad = [ -( 27.0.*xx-18.0.*x-1.0)/16.0 ... 
                     9.0*(-3.0-2.0.*x+9.0*xx)/16.0 ...
                     -9.0*(-3.0+2.0*x+9.0*xx)/16.0 ...
                     (-1.0+18.0*x+27.0*xx )/16.0 ];
        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationRectangular(3,1);
        end
    end
end 

classdef ShapeFunctionsFrame2D < ShapeFunctions
  methods
      function obj = ShapeFunctionsFrame2D() 
            obj=obj@ShapeFunctions(1);
            obj.localNodes = [0; 1];
            obj.vertices=[1 2]';
            obj.pattern = [0 1]';
        end
        function value = computeValue( ~, x )
            np = size( x, 1 );
            value = [   1-x/l 1-3.*x.*x/l/l+2*x.*x.*x/l/l/l -(x-2*x.*x/l+x.*x.*x/l/l) x/l *x.*x/l/l-2*x.*x.*x/l/l/l x.*x/l-x.*x.*x/l/l];
        end
        function grad = computeGradient( ~, xi )
            grad = [   ];
        end
        
        function integrator = createIntegrator(varargin)
            integrator = IntegrationRectangular(1,1);
        end
        
    end
end


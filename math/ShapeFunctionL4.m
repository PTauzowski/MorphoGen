classdef ShapeFunctionL4 < ShapeFunctions2D
    methods
        function obj = ShapeFunctionL4()
            obj.localNodes = [-1 -1; 1 -1; -1 1; 1 1 ];
            obj.vertices=[1 2 3 4]';
            obj.edges = [1 2; 2 4; 4 3; 3 1 ]';
            obj.contour = [1 2 4 3];
            obj.edgesf = ShapeFunctionL2;
            obj.pattern = [ 0 1 0 1; 0 0 1 1 ]';
        end
        function value = computeValue( ~, xi )
            r = xi(:,1);
            s = xi(:,2);
            np = size( xi, 1 );
            value = [ (1-r(1:np)).*(1-s(1:np))/4 (1+r(1:np)).*(1-s(1:np))/4 (1-r(1:np)).*(1+s(1:np))/4 (1+r(1:np)).*(1+s(1:np))/4 ];
        end
        function grad = computeGradient( ~, xi )
            r  = xi(:,1);
            s  = xi(:,2);
            grad = zeros( 4, 2, size( xi, 1 ) );
            grad(1,1,:) = -(1-s)/4;  grad(2,1,:) =  (1-s)/4;  grad(3,1,:) = -(1+s)/4; grad(4,1,:) = (1+s)/4;  
            grad(1,2,:) = -(1-r)/4;  grad(2,2,:) = -(1+r)/4;  grad(3,2,:) =  (1-r)/4; grad(4,2,:) = (1+r)/4; 
        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationRectangular(1,2);
        end
    end
end 

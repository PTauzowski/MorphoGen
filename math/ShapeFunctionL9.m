classdef ShapeFunctionL9 < ShapeFunctions2D
    methods
        function obj = ShapeFunctionL9() 
            obj.localNodes = [-1 -1; 0 -1; 1 -1; -1 0; 0 0; 1 0; -1 1; 0 1; 1 1];
            obj.vertices=[1 3 7 9]';
            obj.edges = [1 2 3; 3 6 9; 9 8 7; 7 4 1]';
            obj.contour = [1 2 3 6 9 8 7 4 1];
            obj.edgesf = ShapeFunctionL3;
            obj.pattern = [ 0 1 2 0 1 2 0 1 2; 0 0 0 1 1 1 2 2 2 ]';
        end
        function value = computeValue( ~, xi )
            x = xi(:,1);
            y = xi(:,2);
            xx = x .* x;
            yy = y .* y;
            value = [ ((x-1.0).*x.*(y-1.0).*y)/4.0 ...
                      ((1.0-xx).*(y-1.0).*y)/2.0 ...
                        (x.*(x+1.0).*(y-1.0).*y)/4.0 ...
                        ((x-1.0).*x.*(1.0-yy))/2.0 ...
                        (1.0-xx).*(1.0-yy) ...
                        (x.*(x+1.0).*(1.0-yy))/2.0 ...
                        ((x-1.0).*x.*y.*(y+1.0))/4.0 ...
                        ((1.0-xx).*y.*(y+1.0))/2.0 ...
                        (x.*(x+1.0).*y.*(y+1.0))/4.0 ];
        end
        function grad = computeGradient( ~, xi )
            x  = xi(:,1);
            y  = xi(:,2);
            grad = zeros( 9, 2, size( xi, 1 ) );
                     
            xx = x.*x;
            yy = y.*y;
            xm1 = x-1.0;
            ym1 = y-1.0;
            xp1 = x+1.0;
            yp1 = y+1.0;

            grad(1,1,:) = (x.*(ym1).*y)/4.0+((xm1).*(ym1).*y)/4.0;
            grad(2,1,:) = -x.*(ym1).*y;
            grad(3,1,:) = ((xp1).*(ym1).*y)/4.0+(x.*(ym1).*y)/4.0;
            grad(4,1,:) = (x.*(1.0-yy))/2.0+((xm1).*(1.0-yy))/2.0;
            grad(5,1,:) = -2.0*x.*(1.0-yy);
            grad(6,1,:) = ((xp1).*(1.0-yy))/2.0+(x.*(1.0-yy))/2.0;
            grad(7,1,:) = (x.*y.*(yp1))/4.0+((xm1).*y.*(yp1))/4.0;
            grad(8,1,:) = -x.*y.*(1.0+y);
            grad(9,1,:) = ((xp1).*y.*(yp1))/4.0+(x.*y.*(yp1))/4.0;

            grad(1,2,:) = ((xm1).*x.*y)/4.0+((xm1).*x.*(ym1))/4.0;
            grad(2,2,:) = ((1.0-xx).*y)/2.0+((1.0-xx).*(ym1))/2.0;
            grad(3,2,:) = (x.*(xp1).*y)/4.0+(x.*(xp1).*(ym1))/4.0;
            grad(4,2,:) = -(xm1).*x.*y;
            grad(5,2,:) = -2.0*(1.0-xx).*y;
            grad(6,2,:) = -x.*(1.0+x).*y;
            grad(7,2,:) = ((xm1).*x.*(yp1))/4.0+((xm1).*x.*y)/4.0;
            grad(8,2,:) = ((1.0-xx).*(yp1))/2.0+((1.0-xx).*y)/2.0;
            grad(9,2,:) = (x.*(xp1).*(yp1))/4.0+(x.*(xp1).*y)/4.0;

        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationRectangular(2,2);
        end
    end
end 

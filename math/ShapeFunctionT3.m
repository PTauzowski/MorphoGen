classdef ShapeFunctionT3 < ShapeFunctions2D
    methods
        function obj = ShapeFunctionT3() 
            obj.localNodes = [ 1 0;  0 1;  0 0 ];
            obj.vertices=[1 2 3]';
            obj.edges = [1 2; 2 3; 3 1 ]';
            obj.contour = [1 2 3];
            obj.edgesf = ShapeFunctionL2;
            obj.pattern = [  ]';
        end
        function value = computeValue( ~, xi )
            value = [ xi(:,1) xi(:,2)  1-xi(:,1)-xi(:,2) ];
        end
        function grad = computeGradient( ~, xi )
            grad = zeros( 3, 2, size( xi, 1 ) );
            grad(1,1,:) = 1;  grad(1,2,:) = 0; 
            grad(2,1,:) = 0;  grad(2,2,:) = 1;  
            grad(3,1,:) = -1;  grad(3,2,:) = -1;  
        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationTriangular(1,2);
        end
        function N = getRecoveryMatrix(obj)
            integrator = obj.createIntegrator();
            N = 3*obj.computeValue( integrator.points )';
        end
    end
end 

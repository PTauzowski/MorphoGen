classdef ShapeFunctionT4 < ShapeFunctions3D
    methods
        function obj = ShapeFunctionT4() 
            obj.localNodes = [ 1 0 0;  0 1 0;  0 0 1;  0 0 0 ];
            obj.vertices=[1 2 3 4]';
            obj.edges = [1 2; 2 4; 4 1; 2 3; 3 1; 3 4 ]';
            obj.faces = [1 2 4; 2 3 4; 3 1 4; 2 1 3]';
            obj.fcontours = [1 2 4; 2 3 4; 3 1 4; 2 1 3]';
            obj.edgesf = ShapeFunctionL2;
            obj.facesf = ShapeFunctionT3;
            obj.pattern = [  ]';
        end
        function value = computeValue( ~, xi )
            value = [ xi(:,1) xi(:,2) xi(:,2) 1-xi(:,1)-xi(:,2)-xi(:,3) ];
        end
        function grad = computeGradient( ~, xi )
            grad = zeros( 4, 3, size( xi, 1 ) );
            grad(1,1,:) = 1;  grad(1,2,:) = 0; grad(1,3,:) = 0; 
            grad(2,1,:) = 0;  grad(2,2,:) = 1; grad(2,3,:) = 0;  
            grad(3,1,:) = 0;  grad(3,2,:) = 0;  grad(3,3,:) = 1; 
            grad(4,1,:) = -1;  grad(4,2,:) = -1;  grad(4,3,:) = -1; 
        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationTriangular(1,3);
        end
        function N = getRecoveryMatrix(obj)
            integrator = obj.createIntegrator();
            N = 3*obj.computeValue( integrator.points )';
        end
    end
end 

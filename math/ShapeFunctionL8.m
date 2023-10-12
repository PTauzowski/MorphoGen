classdef ShapeFunctionL8 < ShapeFunctions3D
    methods
        function obj = ShapeFunctionL8() 
            obj.localNodes = [-1 -1 -1;  1 -1 -1; -1 1 -1; 1 1 -1; ...                            
                              -1 -1  1;  1 -1  1; -1 1  1;  1 1  1];
            obj.vertices = [ 1 2 3 4 5 6 7 8 ]';
            obj.edges    = [ 1 2; 2 4; 4 3; 3 1; 1 5; 2 6; 4 8; 3 7; 5 6; 6 8; 8 7; 7 5 ]';
            obj.edgesf = ShapeFunctionL2;
            obj.faces    = [ 1 2 5 6;  2 4 6 8;  4 3 8 7; 3 1 7 5; 4 3 2 1; 5 6 7 8 ]';
            obj.fcontours    = [ 1 2 6 5;  2 4 8 6; 4 3 7 8; 3 1 5 7; 4 3 1 2; 5 6 8 7 ]';
            obj.facesf = ShapeFunctionL4;
            obj.pattern  = [ 0 1 0 1 0 1 0 1; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1 ]';
        end
        function value = computeValue( ~, xi )
            x = xi(:,1);
            y = xi(:,2);
            z = xi(:,3);
            value = [ ((1-x).*(1-y).*(1-z))/8  ((x+1).*(1-y).*(1-z))/8  ((1-x).*(y+1).*(1-z))/8 ((x+1).*(y+1).*(1-z))/8 ((1-x).*(1-y).*(z+1))/8  ((x+1).*(1-y).*(z+1))/8  ((1-x).*(y+1).*(z+1))/8 ((x+1).*(y+1).*(z+1))/8 ];
        end
        function grad = computeGradient( ~, xi )
             x  = xi(:,1);
             y  = xi(:,2);
             z  = xi(:,3);
             grad = zeros( 8, 3, size( xi, 1 ) );

             grad(1,1,:) = -((1-y).*(1-z))/8;
             grad(2,1,:) =  ((1-y).*(1-z))/8;
             grad(3,1,:) = -((y+1).*(1-z))/8;
             grad(4,1,:) =  ((y+1).*(1-z))/8;
             grad(5,1,:) = -((1-y).*(z+1))/8;
             grad(6,1,:) =  ((1-y).*(z+1))/8;
             grad(7,1,:) = -((y+1).*(z+1))/8;
             grad(8,1,:) =  ((y+1).*(z+1))/8;

             grad(1,2,:) = -((1-x).*(1-z))/8;
             grad(2,2,:) = -((x+1).*(1-z))/8;
             grad(3,2,:) =  ((1-x).*(1-z))/8;
             grad(4,2,:) =  ((x+1).*(1-z))/8;
             grad(5,2,:) = -((1-x).*(z+1))/8;
             grad(6,2,:) = -((x+1).*(z+1))/8;
             grad(7,2,:) =  ((1-x).*(z+1))/8;
             grad(8,2,:) =  ((x+1).*(z+1))/8;

             grad(1,3,:) = -((1-x).*(1-y))/8;
             grad(2,3,:) = -((x+1).*(1-y))/8;
             grad(3,3,:) = -((1-x).*(y+1))/8;
             grad(4,3,:) = -((x+1).*(y+1))/8;
             grad(5,3,:) =  ((1-x).*(1-y))/8;
             grad(6,3,:) =  ((x+1).*(1-y))/8;
             grad(7,3,:) =  ((1-x).*(y+1))/8;
             grad(8,3,:) =  ((x+1).*(y+1))/8;
        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationRectangular(1,3);
        end
    end
end 

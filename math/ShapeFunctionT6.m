classdef ShapeFunctionT6 < ShapeFunctions2D
    methods
        function obj = ShapeFunctionT6() 
            obj.localNodes = [1 0; 0 0.5; 0 1; 0.5 0; 0 0; 0 0];
            obj.vertices=[1 2 3]' ;
            obj.edges = [1 4 2; 2 5 3; 3 6 1 ]';
            obj.contour = [1 4 2 5 3 6];
            obj.edgesf = ShapeFunctionL3;
            obj.pattern = [  ]';
        end
        function value = computeValue( ~, xi )
            L1 = xi(:,1);
            L2 = xi(:,2);
            L3 = 1-L1-L2;
            value = [ (2*L1-1).*L1 ...
                      4*L1.*L2  ...
                      (2*L2-1).*L2  ...
                      4*L2.*L3 ...
                      (2*L3-1).*L3 ...
                      4*L3.*L1 ];
        end
        function grad = computeGradient( ~, xi )
            L1 = xi(:,1);
            L2 = xi(:,2);
            L3 = 1-L1-L2;
            grad = zeros( 3, 2, size( xi, 1 ) );
            
            grad(1,1,:) = 4*L1;  
            grad(2,1,:) = 4*L2;  
            grad(3,1,:) = 0; 
            grad(4,1,:) = -4*L2; 
            grad(5,1,:) = 1-4*L3; 
            grad(6,1,:) = 4*(1-2*L1-4*L2); 
            
            grad(1,1,:) = 0;  
            grad(2,1,:) = 4*L1;  
            grad(3,1,:) = 4*L2; 
            grad(4,1,:) = 4*(1-2*L1-4*L2); 
            grad(5,1,:) = 1-4*L3; 
            grad(6,1,:) = -4*L1; 
%             
        end
        function integrator = createIntegrator(varargin)
            integrator = IntegrationTriangular(2,2);
        end
    end
end 

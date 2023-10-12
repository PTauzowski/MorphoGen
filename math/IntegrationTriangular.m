classdef IntegrationTriangular < Integration
        
    methods
        function obj = IntegrationTriangular(g,d)
            obj.dg = g;
            obj.dim = d;
            
            switch g
                case 0
                    points1D  = 0.0;
                    weights1D = 2.0;
                    
                case 1
                     switch d
                        case 2
                            obj.points = [ 1/3 1/3 ];
                            obj.weights = 1/2;
                        case 3
                            obj.points = [ 1/4 1/4 1/4 ];
                            obj.weights = 1/2;
                         otherwise
                            warning('Unsupported space dimension for linear triangles: ' + num2str( d ) );
                      end
                     
                case 2    
                    switch d
                        case 2
                            obj.points = [ 1/2 1/2; 1/2 0; 0 1/2 ];
                            obj.weights = [ 1/6 1/6 1/6 ];
                        case 3
                            a = 0.58541020;
                            b = 0.13819660;
                            obj.points = [ a b b; b a b; b b a; b b b ];
                            obj.weights = [ 1/8 1/8 1/8 1/8 ];
                         otherwise
                            warning('Unsupported space dimension for quadratic triangles: ' + num2str( d ) );
                      end
                    
                case 3
                    switch d
                        case 2
                            obj.points = [ 1/3 1/3 1/3; 0.6 0.2 0.2; 0.2 0.6 0.2; 0.2 0.2 0.6 ];
                            obj.weights = 0.5*[ -27/48 25/48 25/48 25/48 ];
                        case 3
                            obj.points = [ 1/4 1/4 1/4;  1/2 1/6 1/6;  1/6 1/2 1/6;  1/6 1/6 1/2;  1/6 1/6 1/6];
                            obj.weights = 0.5*[ -4/5 9/20 9/20 9/20 9/20 ];
                         otherwise
                            warning('Unsupported space dimension for cubic triangles: ' + num2str( d ) );
                      end
                otherwise
                    warning('Unsupported triangular quadrature dimension: ' + num2str( g ) );
            end
            
        end
        function vol = integrateVectorFunction( obj, fn ) 
            [dets,~,~] = fn.Jdet( obj.points );
            vol = sum( dets(:) .* obj.weights(:) .* fn.computeValue( obj.points ), 2 );
        end
        function vol = integrateMatrixFunction( obj, fn ) 
            [dets,~,~] = fn.Jdet( obj.points );
            vol = sum( dets(:) .* obj.weights(:) .* fn.computeValue( obj.points ), 3 );
        end
        function vol = integrateScalarFunction( obj, fn ) 
            dets = fn.Jdet( obj.points );
            vol = sum( dets(:) .* obj.weights(:) .* fn.computeValue( obj.points ) );
        end
    end
end


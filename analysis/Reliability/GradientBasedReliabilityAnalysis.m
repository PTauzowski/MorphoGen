classdef GradientBasedReliabilityAnalysis < ReliabilityAnalysis
    
    
    methods
        function obj = GradientBasedReliabilityAnalysis(g, randVars)
            obj = obj@ReliabilityAnalysis(g,randVars);
        end
        
        function [g, dg] = gradG( obj, u )
                dim = obj.getDim();
                dg  = zeros(dim,1);
                x   = u*0;
                pert = 0.0001;

                x = obj.transformToU( u );
                % gn  = LSF( data, xn );
                g = obj.g.computeValue( x );

                for k=1:dim
                        u1 = u;
                        u1( k )  = u1( k ) + pert;
                        x = obj.transformToU( u1 );
                        g = obj.g.computeValue( x );
                        dg( k ) = ( gk - g ) / pert

                end
            end
        end
        
end


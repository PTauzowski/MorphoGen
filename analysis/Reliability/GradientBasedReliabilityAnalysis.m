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
                        gk = u;
                        gk( k )  = gk( k ) + pert;
                        x = obj.transformToU( gk );
                        g = obj.g.computeValue( x );
                        dg( k ) = ( gk - g ) / pert

                end
            end
        end
        
end


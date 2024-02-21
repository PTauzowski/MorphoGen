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

                for subiter=1:10
                  g = obj.g.computeValue( x );
                  if  corr == false
                         x=0.5*x;
                  else
                    break;
                  end
                end

                for k=1:dim
                        u1 = u;
                        u1( k )  = u1( k ) + pert;
                        x = obj.transformToU( u1 );
                        g = obj.g.computeValue( x );
                        if  corr == false
                            dg( k ) = 1;
                        else
                            dg( k ) = ( gk - g ) / pert;
                        end

                end
            end
        end
        
end


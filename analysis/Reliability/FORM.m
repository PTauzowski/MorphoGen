classdef FORM < GradientBasedReliabilityAnalysis
    properties
        betat;
    end
    
    methods
        function obj = FORM(g,randVars)
            obj = obj@GradientBasedReliabilityAnalysis(g,randVars);
        end
        
        function [ Pf, xn, beta ] = solve(obj,x0)
           pert = 0.0001;

            dim = obj.getDim();
            un  = zeros( 1, dim );
            xn  = ones( 1, dim );
            dg  = zeros( 1, dim );
            x1  = zeros( 1, dim );
            convP = 1.0E100;

            for iter=1:100000
                xn = obj.transformFromU( un );

               [gn, corr]  = obj.g.computeValue( xn );

                for k=1:dim

                    u1 = un;
                    u1( k )  = u1( k ) + pert;
                    xl = obj.transformFromU( u1 );

                    [gk, corr] = obj.g.computeValue( xl );
                    %gk = LSF( data, x1 );
                    if  corr == false
                         dg( k ) = 1;
                    else
                        dg( k ) = ( gk - gn ) / pert;
                    end

                end
                if  norm(dg)<1.0E-20
                     Pf = -1;
                     beta = -1;
                     mpp=dg;
                     return;
                end


                unp  = 1.0 / norm( dg )^2 * ( dg * un' - gn ) * dg;
                du   = unp - un;
                conv = norm( du );

                if  size(find( abs(conv) > 50 ),1 ) || iter > 50 
                     Pf = -1;
                     beta = -1;
                     return;
                end

                un   = unp;

                if conv < 0.001 
                    beta = norm( un );
                    Pf = normcdf( -beta );
                    xn = obj.transformFromU( un );
                    return;
                end
                convP=conv;

            end
        end
    end
end


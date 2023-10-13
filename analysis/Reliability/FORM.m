classdef FORM < GradientBasedReliabilityAnalysis
    properties
        betat;
    end
    
    methods
        function obj = FORM(g,randVars)
            obj = obj@GradientBasedReliabilityAnalysis(g,randVars);
        end
        
        function [ Pf, mpp, beta ] = solve(obj)
            pert = 0.0001;

            dim = obj.getDim();
            un  = zeros( 1, dim );
            xn  = ones( 1, dim );
            dg  = zeros( 1, dim );
            x1  = zeros( 1, dim );
            convP = 1.0E100;

            for iter=1:100000
               xn = obj.transformFromU( un );
               [g, dg]  = obj.g.compute( xn );
                
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
                    mpp = obj.transformFromU( un );
                    return;
                end

            end
        end
    end
end


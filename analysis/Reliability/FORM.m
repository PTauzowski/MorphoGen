classdef FORM < GradientBasedReliabilityAnalysis
    properties
        betat;
    end
    
    methods
        function obj = FORM(randVars,g,transformU)
            obj = obj@GradientBasedReliabilityAnalysis(randVars,g,transformU);
        end
        
        function [ Pf, mpp, beta, success ] = solve(obj)
            dim = obj.getDim();
            u  = zeros( 1, dim );
            mpp=u;
            for iter=1:100
               [g, dg] = obj.computeGu( u );
                
               if  norm(dg)<1.0E-20
                     Pf = -1;
                     beta = -1;
                     mpp=dg;
                     success=false;
                     return;
               end
              
                unext  = 1.0 / norm( dg )^2 * ( dg * u' - g ) * dg;
                du   = u - unext;
                conv = norm( du );

                if  size(find( abs(conv) > 50 ),1 )  
                     Pf = -1;
                     beta = -1;
                     mpp=dg;
                     success=false;
                     return;
                end
                u = unext;
                if conv < 0.0001 
                    beta = norm( u );
                    Pf = normcdf( -beta );
                    mpp = obj.transform.toX( u );
                    success=true;
                    return;
                end
            end
            fprintf('FORM not converged!\n');
            Pf = -1;
            beta = -1;
            mpp=dg;
            success=false;
        end
    end
end


classdef FORM < GradientBasedReliabilityAnalysis
    properties
        betat;
    end
    
    methods
        function obj = FORM(randVars,g,transformU)
            obj = obj@GradientBasedReliabilityAnalysis(randVars,g,transformU);
        end
        
        function results = solve(obj,x0)
            dim = obj.getDim();
            u=obj.transform.toU(x0);
            results.g0 = obj.g.computeValue( x0 );
            mpp = u;
            for iter=1:100
               [g, dg] = obj.computeGu( u );
                
               if  norm(dg)<1.0E-20
                     results.success=false;
                     results.err_msg='FORM error: gradient norm too small';
                     return;
               end
              
                unext  = 1.0 / norm( dg )^2 * ( dg * u' - g ) * dg;
                du   = u - unext;
                conv = norm( du );

                if  size(find( abs(conv) > 500 ),1 ) 
                     results.success=false;
                     results.err_msg='FORM not converged!\n';
                     return;
                end
                u = unext;
                if conv < 0.00001 
                    results.beta = norm( u );
                    results.Pf = normcdf( -results.beta );
                    results.mpp = obj.transform.toX( u );
                    results.gmpp = obj.g.computeValue(results.mpp );
                    results.success=true;
                    return;
                end
            end
            results.success=false;
            results.err_msg='FORM not converged!\n';
        end
    end
end


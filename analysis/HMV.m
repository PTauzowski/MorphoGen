classdef HMV < GradientBasedReliabilityAnalysis
    
     properties
        betat;
    end
    
    methods
        function obj = HMV(randVars,g,transform,betat)
            obj = obj@GradientBasedReliabilityAnalysis(randVars,g,transform);
            obj.betat=betat;
        end
        
        function results = solve(obj,x0)
            dim = obj.getDim();
            u   = obj.transform.toU(x0);
            n   = zeros( dim, 1 );
            n1  = zeros( dim, 1 );
            [g, dg] = computeGu( obj, u ); 
            if  norm(dg)<1.0E-20
                     results.success=false;
                     results.err_msg='HMV error: gradient norm too small';
                     return;
            end
            g0=g;
            for k=1:100
                g1 = g;
                n2 = n1; 
                n1 = n;  
                [g, dg] = computeGu( obj, u );
                n = -dg / norm(dg);
                beta = norm( u );
                dgrel = abs((g-g1)/norm(dg) );
                dgabs = abs((g-g1) );
                if  max( abs(beta-obj.betat), max(dgrel, dgabs) ) < 0.0001  
                    
                     betar = obj.betat*(1+g/(g0-g));
                     ur = obj.betat*n;
                     mpp = obj.transform.toX( ur );
                     g = obj.g.computeValue( mpp );

                     results.g0=g0;
                     results.mpp = obj.transform.toX(u);
                     results.mpp0_pred = obj.transform.toX(n*betar);
                     results.g = obj.g.computeValue( results.mpp );
                     results.n=n;
                     results.beta_pred = betar;
                     results.success=true;
                    return;
                end

                if k<3
                    u = obj.betat * n;
                else
                    s = dot( n1 - n , n - n2 );
                    if s>=0
                        u = obj.betat * n;
                    else
                        ne=n+n1+n2;
                        u = obj.betat * ne / norm(ne);
                    end
                end

            end
            results.success = false;
            results.err_msg = ['HMV error: not convergent after ' num2str(k-1) ' iterations'];
        end

        function printResults(obj, tx, hmv_results)
            fprintf("%s HMV results :",tx);
            if hmv_results.success
                fprintf("G0=%5.7f, G=%5.7f, beta_pred=%5.7f, mpp=",hmv_results.g0,hmv_results.g,hmv_results.beta_pred);
                obj.printPoint(hmv_results.mpp);
            else
                fprintf("NOT succeed! %s",hmv_results.err_msg);
            end
            fprintf("\n");
        end
    end
end


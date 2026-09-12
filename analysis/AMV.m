classdef AMV < GradientBasedReliabilityAnalysis
    
    properties
        betat;
    end
    
    methods
        function obj = AMV(randVars,g,transform,betat)
            obj = obj@GradientBasedReliabilityAnalysis(randVars,g,transform);
            obj.betat=betat;
        end
        
        function results = solve(obj,x0)
            dim = obj.getDim();
            u   = obj.transform.toU(x0);
            [g, dg] = computeGu( obj, u );
            if  norm(dg)<1.0E-20
                     results.n = -1;
                     results.mpp=dg;
                     results.success=false;
                     results.err_msg='AMV error: gradient norm too small';
                     return;
            end
            g0=g;
            n = -dg ./ norm(dg);
            u = obj.betat * n;
            for k=1:100
                g1 = g;
                [g, dg] = computeGu( obj, u );
                n = -dg ./ norm(dg,2);
                u = obj.betat * n;
                beta = norm(u);
                dgrel = abs((g-g1)/norm(dg) );
                dgabs = abs(g-g1);
                conv = max( abs(beta-obj.betat), max(dgrel, dgabs) );
                if  conv < 0.0001 
                     results.n = n;
                     results.mpp = obj.transform.toX( u );
                     results.beta_pred = obj.betat*(1+g/(g0-g));
                     results.g0=g0;
                     results.g = g;
                     results.success=true;
                    return;
                end
            end
            results.success = false;
            results.err_msg = ['AMV error: not convergent after ' num2str(k-1) ' iterations'];
        end

        function printResults(obj, tx, amv_results)
            fprintf("%s AMV results :",tx);
            if amv_results.success
                fprintf("G0=%5.7f, G=%5.7f, beta_pred=%5.7f, mpp=",amv_results.g0,amv_results.g,amv_results.beta_pred);
                obj.printPoint(amv_results.mpp);
            else
                fprintf("NOT succeed! %s",amv_results.err_msg);
            end
            fprintf("\n");
        end
    end
end


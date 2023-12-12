classdef AMV < GradientBasedReliabilityAnalysis
    
    properties
        betat;
    end
    
    methods
        function obj = AMV(randVars,g,transform,betat)
            obj = obj@GradientBasedReliabilityAnalysis(randVars,g,transform);
            obj.betat=betat;
        end
        
        function [ Pf, mpp, n, success ] = solve(obj)
            dim = obj.getDim();
            u   = zeros( 1, dim );
            [g, dg] = computeGu( obj, u );

            if  norm(dg)<1.0E-20
                     Pf = -1;
                     n = -1;
                     mpp=dg;
                     success=false;
                     return;
            end

            g0=g;
            n = -dg ./ norm(dg);
            u = obj.betat * n;
            mpp=u;
            for k=1:100
                g1 = g;
                [g, dg] = computeGu( obj, u );
                n = -dg ./ norm(dg,2);
                u = obj.betat * n;
                %obj.transform.toX( u )
                beta = norm(u);
                dgrel = abs((g-g1)/norm(dg) );
                dgabs = abs(g-g1);

                conv = max( abs(beta-obj.betat), max(dgrel, dgabs) );
                betar = obj.betat*(1+g/(g0-g));
                if  conv < 0.0001 
                     Pf = normcdf( -norm(u) );
                     ur = betar*n;
                     mpp = obj.transform.toX( ur );
                     success=true;
                    return;
                end
            end
            fprintf('AMV not converged!\n');
            Pf=0;
            beta=0.5;
        end
    end
end


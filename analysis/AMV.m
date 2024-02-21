classdef AMV < GradientBasedReliabilityAnalysis
    
    properties
        betat;
    end
    
    methods
        function obj = AMV(g,randVars,betat)
            obj = obj@GradientBasedReliabilityAnalysis(g,randVars);
            obj.betat=betat;
        end
        
        function [ Pf, mpp, beta ] = solve(obj,x0)
            dim = obj.getDim();
            u   = zeros( dim, 1 );
            epsilon=1.0e-03;
            g=0;
            [g, dg] = obj.gradG( u );
            if  norm(dg)<1.0E-20
                     Pf = -1;
                     beta = -1;
                     mpp=dg;
                     return;
            end
            g0=g;
            n = -dg ./ norm(dg);
            u = obj.betat * n;
            for k=1:100000
                g1 = g;
                [g, dg] = obj.gradG( u );
                n = -dg ./ norm(dg,2);
                u = obj.betat * n;
                beta = norm(u);
                dgrel = abs((g-g1)/norm(g) );
                dgabs = abs(g-g1);

                conv = max( abs(beta-obj.betat), max(dgrel, dgabs) );
                betar = obj.betat*(1+g/(g0-g));
                if  conv < epsilon 
                     Pf = normcdf( -norm(u) );
                     ur = betar*n;
                     mpp = obj.transformFromU( ur );
                     %[g, ~]  = obj.g.computerValue( mpp );
                     mpp=n;
                     %fprintf('G = %5.3f\n',g);
                    return;
                end
             end
        end
    end
end


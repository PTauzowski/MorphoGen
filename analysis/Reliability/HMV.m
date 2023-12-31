classdef HMV < GradientBasedReliabilityAnalysis
    
     properties
        betat;
    end
    
    methods
        function obj = HMV(g,randVars,betat)
            obj = obj@GradientBasedReliabilityAnalysis(g,randVars);
            obj.betat=betat;
        end
        
        function [ Pf, mpp, betar ] = solve(obj,x0)
            dim = obj.getDim();
            u   = zeros( dim, 1 );
            n   = zeros( dim, 1 );
            n1  = zeros( dim, 1 );
            dg  = zeros( dim, 1 );
            epsilon=1.0e-03;
            [g, dg] = obj.gradG( u );  
            if  norm(dg)<1.0E-20
                     Pf = -1;
                     betar = -1;
                     mpp=dg;
                     return;
            end
            g0=g;
            for k=1:100000
                g1 = g;
                n2 = n1; 
                n1 = n;  
                [g, dg] = obj.gradG( u );
                n = -dg ./ norm(dg);
                beta = norm( u );
                dgrel = abs((g-g1)/norm(g) );
                dgabs = abs((g-g1) );
                if  max( abs(beta-obj.betat), max(dgrel, dgabs) ) < epsilon 
                     %Pf = normcdf( -norm( u ) );
                     mpp = zeros(dim,1);
                     betar = obj.betat*(1+g/(g0-g));
                     ur = betar*n;
                     mpp = obj.transformFromU( ur );
                     [g, corr]  = obj.g.computeValue( mpp );
                     mpp = n;
                     Pf = normcdf( -beta );
                     % fprintf('iter = %2d\n',k);
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
        end
    end
end


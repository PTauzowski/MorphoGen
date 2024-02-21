classdef CMV < GradientBasedReliabilityAnalysis
    
    properties
        betat;
    end
    
    methods
        function obj = CMV(g,randVars,betat)
            obj = obj@GradientBasedReliabilityAnalysis(g,randVars);
            obj.betat=betat;
        end
       
        function [ Pf, mpp, beta ] = solve(obj,x0)
            dim = obj.getDim();
            u   = zeros( dim, 1 );
            n   = zeros( dim, 1 );
            n1  = zeros( dim, 1 );
            epsilon=1.0e-03;
            [g, dg] = obj.gradG( u );
            if  norm(dg)<1.0E-20
                     Pf = -1;
                     beta = -1;
                     mpp=dg;
                     return;
            end
            for k=1:100000
                g1 = g;
                n2 = n1; 
                n1 = n;  
                [g, dg] = obj.gradG( u );
                %fprintf('G = %5.3f\n',g);
                n = -dg ./ norm(dg);
                beta = norm(u);
                dgrel = abs((g-g1)/norm(g));
                dgabs = abs((g-g1));

                if  max( abs(beta-obj.betat), max(dgrel, dgabs) ) < epsilon 
                     Pf = normcdf(-norm(u));
                     mpp = zeros(dim,1);
                     mpp = obj.transformFromU( u );
                     [g, ~]  = obj.g.computeValue( mpp );
                     %fprintf('G = %5.3f\n',g);
                    return;
                end

                if k<3
                    u = obj.betat * n;
                else
                    s = dot( n1 - n , n - n2 );
                    ne=n+n1+n2;
                    u = obj.betat * ne / norm(ne);
                end

            end
        end
    end
end


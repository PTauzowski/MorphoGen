classdef HMV < GradientBasedReliabilityAnalysis
    
     properties
        betat;
    end
    
    methods
        function obj = HMV(randVars,g,betat)
            obj = obj@GradientBasedReliabilityAnalysis(randVars,g);
            obj.betat=betat;
        end
        
        function [ Pf, mpp, betar ] = solve(obj,x0)
            dim = obj.getDim();
            J   = zeros(dim);
            u   = zeros( 1, dim );
            n   = zeros( 1, dim );
            n1  = zeros( 1, dim );
            dg  = zeros( 1, dim );
            epsilon=1.0e-03;
            x = obj.transformFromU( u );
            [g, dg] = obj.g.compute( x );
            for k=1:size(obj.randVars,2)
                        J(k,k)=normpdf(u(k))/pdf(obj.randVars{k}.pd,x(k));
            end
            dg=dg*J;
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
                x = obj.transformFromU( u );
                [g, dg] = obj.g.compute( x );
                for j=1:size(obj.randVars,2)
                        J(j,j)=normpdf(u(j))/pdf(obj.randVars{k}.pd,x(j));
                end
                dg=dg*J;
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
                     g = obj.g.computeValue( mpp );
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


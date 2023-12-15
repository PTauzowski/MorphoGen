classdef HMV < GradientBasedReliabilityAnalysis
    
     properties
        betat;
    end
    
    methods
        function obj = HMV(randVars,g,transform,betat)
            obj = obj@GradientBasedReliabilityAnalysis(randVars,g,transform);
            obj.betat=betat;
        end
        
        function results = solve(obj)
            dim = obj.getDim();
            u   = zeros( dim, 1 );
            n   = zeros( dim, 1 );
            n1  = zeros( dim, 1 );
            [g, dg] = computeGu( obj, u ); 
            if  norm(dg)<1.0E-20
                     results.Pf = -1;
                     results.n = -1;
                     results.mpp=dg;
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
                     %Pf = normcdf( -norm( u ) );
                     %mpp = zeros(dim,1);
                     %betar = norm(u);
                     betar = obj.betat*(1+g/(g0-g));
                     ur = betar*n;
                     mpp = obj.transform.toX( u );
                     g = obj.g.computeValue( mpp );
                     mpp = u;
                     Pf = normcdf( -beta );
                     % fprintf('iter = %2d\n',k);
                     results.mpp = obj.transform.toX( ur );
                     results.g = obj.g.computeValue( mpp );
                     results.n=n;
                     results.beta_pred = obj.betat*(1+g/(g0-g));
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
    end
end


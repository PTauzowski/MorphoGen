classdef GradientBasedReliabilityAnalysis < ReliabilityAnalysis  
    
    properties
        transform;
    end

    methods
        function obj = GradientBasedReliabilityAnalysis(randVars, g, transform)
            obj = obj@ReliabilityAnalysis(randVars,g);
            obj.transform = transform;
            g.setPerturbation( transform.createXPerturbation(0.0001) );
        end
        
        function [g, gradU] = computeGu( obj, u  )
            pert=0.0001;
            x = obj.transform.toX ( u );
            x1 = obj.transform.toX ( u+pert );
            dgU=(x1-x)/pert;
            [g, dgX]  = obj.g.compute( x );
            gradU = obj.transform.gradientToU(x,u,dgX);
        end
    end
        
end


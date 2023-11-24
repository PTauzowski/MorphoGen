classdef GradientBasedReliabilityAnalysis < ReliabilityAnalysis  
    
    properties
        transform;
    end

    methods
        function obj = GradientBasedReliabilityAnalysis(randVars, g, transform)
            obj = obj@ReliabilityAnalysis(randVars,g);
            obj.transform = transform;
            g.setPerturbation( transform.createXPerturbation(0.00001) );
        end
        
      

        function [g, gradU] = computeGu( obj, u  )
            x = obj.transform.toX ( u );
            [g, dgX]  = obj.g.compute( x );
            gradU = obj.transform.gradientToU(x,u,dgX);
        end
    end
        
end


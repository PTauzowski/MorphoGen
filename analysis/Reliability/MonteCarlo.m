classdef MonteCarlo < ReliabilityAnalysis
        
    properties
        nsamples;
        x;
        r;
    end
    
    methods
        function obj = MonteCarlo(g,randVars,nsamples)
            obj = obj@ReliabilityAnalysis(g,randVars);
            obj.nsamples=nsamples;
        end
        
        function [ Pf, r ] = solve(obj)
            obj.x=obj.generateRandomSapmles(obj.nsamples);
            r =  obj.g.computeValue( obj.x );
            obj.r = r;
            Pf = max(size(find(obj.r<0)))/obj.nsamples;
        end
    end
end


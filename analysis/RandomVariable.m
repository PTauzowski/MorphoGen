classdef RandomVariable < handle
    
    properties
        pd,name,mean,sd;
    end
        
    methods
        function obj = RandomVariable(name,m,s)
            obj.pd = makedist(name,m,s);
            obj.mean=m;
            obj.sd=s;
            obj.name=name;
        end
        
        function x = randomize(obj,nsamples)
            x = random(obj.pd,nsamples,1);
        end

        function setMean(obj,m)
            obj.pd = makedist(obj.name,m,obj.sd);
            obj.mean=m;
        end

        function u = toU(obj,x)
            u = norminv( cdf( obj.pd, x ) );
        end

        function x = fromU(obj,u)
               %x = norminv( cdf(obj.pd, u), obj.pd.mean, obj.pd.sigma );
               % if u < 0 
               %      x = norminv( normcdf( max(-38, u) ), obj.mean, obj.sd );
               % else
               %      x = abs(norminv( normcdf( max(-38, -u) ))) * obj.sd + obj.mean;
               % end
               x=icdf(obj.pd, normcdf(u));
               
        end
             
    end
end


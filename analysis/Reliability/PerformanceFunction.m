classdef (Abstract) PerformanceFunction < Function

    properties
        randVars,J;
    end
   
    methods(Abstract)
        function obj = PerformanceFunction(randVars)
            obj=obj@Function(size(randVars,2),0.00001)
            obj.randVars = randvars;
            for k=1:size(obj.randVars,2)
                        obj.eps(k)=obj.eps(k)*obj.randVars{k}.pd.sigma;
            end
            obj.J=zeros(obj.dim,obj.dim);
        end

        function uG = transformGradToU(obj, x, u, dg)
            
        end

        function [g, dgX] = computeU(obj, u)
            [g dgX] = comput
        end

       
    end

end


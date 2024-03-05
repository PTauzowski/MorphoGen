classdef StandardSpaceTransformation < handle
  
    properties
        J,randVars;
    end

    methods (Abstract)
        u = toU(x);
        x = toX(u);
        udG = gradientToU(u,x,dg);
        epsX = createXPerturbation(obj,epsU);
    end
    
    methods
        function obj = StandardSpaceTransformation(randVars)
            obj.randVars = randVars;
        end
    end
end


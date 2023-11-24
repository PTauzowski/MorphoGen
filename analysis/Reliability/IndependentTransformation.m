classdef IndependentTransformation < StandardSpaceTransformation
    
    methods
        function obj = IndependentTransformation(randVars)
            obj = obj@StandardSpaceTransformation(randVars);
            dim=size(obj.randVars, 2);
            obj.J=zeros(dim,dim);
        end
       

        function u = toU(obj,x)
            dim = size(obj.randVars, 2);
            u = zeros( 1, dim );
            for k=1:dim
               u(k)=obj.randVars{k}.toU(x(k)); 
            end      
        end

        function x = toX(obj,u)
            dim = size(obj.randVars, 2);
            x = zeros( 1, dim );
            for k=1:dim
               x(k)=obj.randVars{k}.fromU(u(k)); 
            end
        end

        function udG = gradientToU(obj,x,u,dg)
            for k=1:size(obj.randVars,2)
                        obj.J(k,k)=normpdf(u(k))/pdf(obj.randVars{k}.pd,x(k));
            end
            udG=dg*obj.J;
        end

        function epsX = createXPerturbation(obj, epsU)
            dim = size(obj.randVars, 2);
            epsX=zeros(1,dim);
            for k=1:size(obj.randVars,2)
                        epsX(k)=epsU*obj.randVars{k}.pd.sigma;
            end
        end

    end
end


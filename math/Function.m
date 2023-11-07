classdef (Abstract) Function < handle
   
    properties
        dim, eps, dx;
    end
    methods (Abstract)
        computeValue( x );
    end
    methods
        function obj = Function(dim,eps)
            obj.dim=dim;
            obj.setPerturbation(eps);            
        end
        function [value, grad] = compute( obj, x )
                value = computeValue(obj, x);
                grad  = computeGradient(obj,x);
        end
        function setPerturbation( obj, eps )
            if size(eps,2) == 1
                obj.eps=zeros(obj.dim,1);
                obj.eps=eps;
            elseif size(eps,2)==obj.dim
                %obj.eps=diag(eps);
                obj.eps=eps';
            else
                error(['Incompatibile dimensions. Function dimension:' num2str(obj.dim) ' perturbation dimmension :' num2str(size(eps,2))]);
            end
            obj.dx=diag(obj.eps);
        end
        function grad = computeGradient( obj, x )
            n=size(x,1);
            y=obj.computeValue(x);
            grad=reshape((obj.computeValue(repelem(x,obj.dim,1)+repmat(obj.dx,n,1))-repmat(y,n,1))./repmat(obj.eps,n,1),n,obj.dim);
        end
    end
end


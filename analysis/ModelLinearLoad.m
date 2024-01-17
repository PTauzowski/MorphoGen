classdef ModelLinearLoad < ModelLinear
    properties
        P0fem, u0fem, ures, E, nu;
    end
    
    methods
      
        function setOneX(obj)
            obj.setX(ones(size(obj.x)));
        end

        function obj = setX(obj,x)
            obj.x=x;
            obj.u0fem = obj.analysis.solveWeighted(obj.x);
            dim=size(obj.u0fem,2);
            obj.ures=zeros(1,dim);
            for k=1:dim
                unodal=obj.analysis.fromFEMVector(obj.u0fem(:,k));
                obj.ures(k)=unodal(obj.result_node,k);
            end
        end

        function u = computeLinearDisplacement(obj,rvr)
            u=rvr*obj.ures';
        end
    
    end
end


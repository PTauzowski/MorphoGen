
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
            if size(obj.P0fem,1)==0
                obj.P0fem=obj.analysis.Pfem;
            else
                obj.analysis.Pfem=obj.P0fem;
                obj.analysis.Pnodal(:)=0;
            end
            obj.u0fem = obj.analysis.solveWeighted(obj.x);
            dim=size(obj.u0fem,2);
            obj.ures=zeros(1,dim);
            for k=1:dim
                unodal=obj.analysis.fromFEMVector(obj.u0fem(:,k));
                obj.ures(k)=unodal(obj.result_node,2);
            end
        end

        function u = computeLinearDisplacement(obj,rvr)
            u=rvr*obj.ures';
        end

        function sHM = computeLinearHMstress(obj,pressure)
            np=size(pressure,1);
            sHM=zeros(np,1);
            for k=1:np
                obj.setDisplacement(pressure(k,:));
                obj.analysis.computeElementResults(obj.x);
                sHM(k)=obj.fe.results.nodal.all(obj.result_node,obj.result_number);
            end
        end 

        function setDisplacement(obj,rvr)
            dim=size(rvr,2);
            obj.analysis.qfem=zeros(obj.analysis.getTaskDim(),1);
            for k=1:dim
                obj.analysis.qfem=obj.analysis.qfem+rvr(k)*obj.u0fem(:,k);
            end
            obj.analysis.qnodal=obj.analysis.fromFEMVector(obj.analysis.qfem);
        end
    
    end
end


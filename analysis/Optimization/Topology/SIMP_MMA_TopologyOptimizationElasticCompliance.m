classdef SIMP_MMA_TopologyOptimizationElasticCompliance < SIMP_MMA_TopologyOptimization
    
    properties
        VolConstr;
    end
    
    methods

        function obj = SIMP_MMA_TopologyOptimizationElasticCompliance(Rmin,problem,penal,VolConstr,is_const)
            obj = obj@SIMP_MMA_TopologyOptimization(1,Rmin,problem,penal,is_const);
            tne=obj.FEproblem.getTotalElemsNumber();
            obj.gradConstrValues = zeros(1,1);
            obj.V0=tne*VolConstr;
            obj.VolConstr=VolConstr;
            obj.x(1:obj.totalFENumber,1)=1;
            obj.is_const=is_const;
        end

        function computeObjectiveFunctonWithGradient( obj, x )
            tne=obj.FEproblem.getTotalElemsNumber();
            c=0;
            dc = zeros(tne,1);
            x_ones = ones(tne,1);
            ind=1;
            for i=1:size(obj.FEproblem.felems,1)
               if obj.is_const
                   fsName='computeStifnessMatrixConst';
               else
                   fsName='computeStifnessMatrix';
               end
               nelems = size(obj.FEproblem.felems{i}.elems,1);
               nnodes = size(obj.FEproblem.felems{i}.elems,2);
               ndofs = size( obj.FEproblem.felems{i}.ndofs,2);
               dim = nnodes*ndofs;
               K = reshape(obj.FEproblem.felems{i}.(fsName)(obj.FEproblem.mesh.nodes,x_ones),dim,dim,nelems);
               obj.qnodal = obj.FEproblem.solveWeighted((obj.x).^obj.penal);
               obj.FEproblem.computeElementResults(obj.qnodal,obj.x);
               qelems = obj.FEproblem.felems{i}.createElemSolutionVectors(obj.qnodal);
               for j=1:nelems
                    c = c + obj.x(ind)^obj.penal*qelems(:,j)'*K(:,:,j)*qelems(:,j);
                    dc(ind) = -obj.penal*obj.x(ind)^(obj.penal-1)*qelems(:,j)'*K(:,:,j)*qelems(:,j);
                    ind=ind+1;
               end
            end
            obj.FobjValue = c;
            obj.gradFobjValue = dc;
        end

        function computeConstraintsAndGradient( obj, x)
            obj.constrValues = sum(x(:))/obj.V0-1;
            obj.gradConstrValues(1:obj.totalFENumber,1) = 1.0/obj.V0;
        end

        function printIterationInfo(obj) 
            fprintf('%5i ',obj.iteration );
            fprintf('Vrel=%5.2f ',round(obj.VolConstr*sum( obj.x )/obj.V0*1000)/10);
            fprintf('constr=%6.4f ',obj.constrValues);
            fprintf('change=%6.4f ',obj.change);
            fprintf('\n');
        end
        
    end
end


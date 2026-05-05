classdef LinearElasticityWeighted < FEAnalysis
   
   properties
        isConst;
        cachedWeightedX;
        cachedWeightedKvals;
        cachedWeightedFunction;
   end
   
   methods       
       function obj = LinearElasticityWeighted(felems, mesh, isConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.isConst=isConst;
            obj.rotations=[];
       end
       function qfem = solveWeighted(obj, x, retainStiffness)
           if nargin < 3
               retainStiffness = false;
           end
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           stiffnessFunction = obj.weightedStiffnessFunction();
           Kvals = obj.globalMatrixAggregationWeighted(stiffnessFunction, x);
           if retainStiffness
                obj.cachedWeightedX = x(:);
                obj.cachedWeightedKvals = Kvals;
                obj.cachedWeightedFunction = stiffnessFunction;
           else
                obj.clearWeightedStiffnessCache();
           end
           obj.qfem = solver.solve(Kvals, obj.Pfem);
           qfem=obj.qfem;
           obj.qnodal=obj.fromFEMVector(qfem(:,1));
       end

       function lambda = solveAdjointWithLoad(obj, xPenal, P_adj_fem)
           % Solve K(xPenal)*lambda = P_adj_fem. K is symmetric so the adjoint
           % system is identical to the forward system with a different RHS.
           % P_adj_fem: nDof x nAdjoints — MATLAB sparse \ handles multiple RHS.
           [I, J, ~] = obj.globalMatrixIndices();
           if size(obj.rotations, 1) == 0
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports), obj.rotations);
           end
           stiffnessFunction = obj.weightedStiffnessFunction();
           if obj.hasCachedWeightedStiffness(xPenal, stiffnessFunction)
               Kvals = obj.cachedWeightedKvals;
           else
               Kvals = obj.globalMatrixAggregationWeighted(stiffnessFunction, xPenal);
           end
           lambda = solver.solve(Kvals, P_adj_fem);
           obj.clearWeightedStiffnessCache();
       end

       function stiffnessFunction = weightedStiffnessFunction(obj)
           if obj.isConst
               stiffnessFunction = 'computeStifnessMatrixConst';
           else
               stiffnessFunction = 'computeStifnessMatrix';
           end
       end

       function tf = hasCachedWeightedStiffness(obj, x, stiffnessFunction)
           tf = ~isempty(obj.cachedWeightedKvals) && ...
               isequal(obj.cachedWeightedFunction, stiffnessFunction) && ...
               numel(obj.cachedWeightedX) == numel(x) && ...
               isequal(obj.cachedWeightedX, x(:));
       end

       function clearWeightedStiffnessCache(obj)
           obj.cachedWeightedX = [];
           obj.cachedWeightedKvals = [];
           obj.cachedWeightedFunction = '';
       end

       function obj = saveMatrices(obj,filename)
          [I,J,~] = obj.globalMatrixIndices();
          if obj.isConst
                K = sparse(I,J,obj.globalMatrixAggregationWeighted('computeStifnessMatrixConst',1));
          else
                K = sparse(I,J,obj.globalMatrixAggregationWeighted('computeStifnessMatrix',1));
          end
           supports = obj.supports;
           obj.prepareRHSVectors();
           P = obj.Pfem;
           nodes = obj.mesh.nodes;
           elems = obj.mesh.elems;
           save(filename, "nodes","elems","K","P","supports");
       end
   end
end

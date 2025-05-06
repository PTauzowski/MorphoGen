classdef LinearStability < FEAnalysis

    properties
        lambdas, qforms, freedofs, fixeddofs;
    end

   methods
       function obj = LinearStability(felems, mesh)
            obj= obj@FEAnalysis( felems, mesh );
            obj.rotations=[];
       end
       
       function solve(obj, num_eigenvalues)
           [I,J,~,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
           if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
               R=0;
           end
           obj.freedofs=solver.freedofs;
           obj.fixeddofs=solver.supdofs;
           K = obj.assemblyGlobalMatrix('computeStifnessMatrix',1,false);
           obj.qfem = solver.solve(K, obj.Pfem );
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults();
           Kg = obj.assemblyGlobalMatrix('computeGeometricStifnessMatrix',1,false);
           [obj.qforms, obj.lambdas]=solver.solveEigenproblem(K,Kg,num_eigenvalues);
       end

       function solveWeighted (obj, x, num_eigenvalues )
           [I,J,~,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
           if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
               R=0;
           end
           obj.freedofs=solver.freedofs;
           K = obj.globalMatrixAggregationWeighted('computeStifnessMatrix',x);
           obj.qfem = solver.solve(K, obj.Pfem );
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults(x);
           Kg = obj.globalMatrixAggregationWeighted('computeGeometricStifnessMatrix',x);
           [obj.qforms,obj.lambdas]=solver.solveEigenproblem(K,Kg,num_eigenvalues);
       end

       function setForm(obj,i)
           obj.qfem =0*obj.Pfem;
           obj.qfem(obj.freedofs)= obj.qforms(:,i);
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults();
       end
   end
end


classdef LinearElasticity < FEAnalysis

   methods
       function obj = LinearElasticity(felems, mesh)
            obj= obj@FEAnalysis( felems, mesh );
            obj.rotations=[];
       end
       
       function [qn, K] = solve(obj)
           [I,J,~,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
           if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
               R=0;
           end
           K = obj.globalMatrixAggregation('computeStifnessMatrix');
           obj.qfem = solver.solve(K, obj.Pfem );
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults(obj.qnodal);
           qn=obj.qnodal;
       end
   end
end


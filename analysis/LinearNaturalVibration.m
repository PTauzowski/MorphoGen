classdef LinearNaturalVibration < FEAnalysis

    properties
        lambdas, qforms, freedofs, fixeddofs;
    end

   methods
       function obj = LinearNaturalVibration(felems, mesh)
            obj= obj@FEAnalysis( felems, mesh );
            obj.rotations=[];
       end

       function K = globalMatrixAggregationWeighted(obj, fname, x)
            K = [];
            ei = getElemIndices(obj);
            for k=1:size(obj.felems,2)
                if ismethod(obj.felems{k},fname)
                    K = [ K; obj.felems{k}.(fname)(obj.mesh.nodes,x(ei{k})) ];
                else
                    error("Class " + class(obj.felems{k}) + " or its predecessors not implements function :"+fname);
                end
            end
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
           K = obj.globalMatrixAggregation('computeStifnessMatrix');
           %obj.qfem = solver.solve(K, obj.Pfem );
           %obj.qnodal = obj.fromFEMVector( obj.qfem );
           %obj.computeElementResults();
           M = obj.globalMatrixAggregation('computeMassMatrix');
           [obj.qforms, obj.lambdas]=solver.solveEigenproblem(K,M,num_eigenvalues);
       end

       function solveWeighted(obj, x, num_eigenvalues)
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
           %obj.qfem = solver.solve(K, obj.Pfem );
           %obj.qnodal = obj.fromFEMVector( obj.qfem );
           %obj.computeElementResults(x);
           M = obj.globalMatrixAggregationWeighted('computeMassMatrix',x);
           [obj.qforms,obj.lambdas]=solver.solveEigenproblem(K,M,num_eigenvalues);
       end

       function setForm(obj,i)
           obj.qfem =0*obj.Pfem;
           obj.qfem(obj.freedofs)= obj.qforms(:,i);
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults();
       end
   end
end


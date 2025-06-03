classdef ElasticHarmonicVibrations < FEAnalysis
   
   properties
        isConst,lambda, omegas, p;
   end
   
   methods       
       
       function obj = ElasticHarmonicVibrations(felems, mesh, p, isConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.p = p;
            obj.isConst=isConst;
            obj.rotations=[];
       end

       function [K, M] = computeMatrices(obj,x)
           K  = obj.createSparseMatrix( obj.assemblyGlobalMatrix('computeStifnessMatrix', x, obj.isConst) );
           M  = obj.createSparseMatrix( obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isConst) );
       end
       
       % Harmonic Vibrations
       % function qfem = solve(obj, x)
       %     [I,J,~] = obj.globalMatrixIndices();
       %     obj.prepareRHSVectors();
       %      if size(obj.rotations,1)== 0 
       %         solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
       %     else
       %         solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
       %     end
       %     K  = obj.assemblyGlobalMatrix('computeStifnessMatrix', x, obj.isConst);
       %     M  = obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isConst);
       %     [mode, lambdas] = solver.solveEigenproblem(K,M,10);
       %     lambdas = diag(lambdas(1:5,1:5));
       %     omegas=sqrt(lambdas);
       %     %obj.qfem = solver.solve(K-(lambdas(3)+lambdas(4))/2*M, obj.Pfem);
       %     obj.qfem = solver.solve(K-0.5*(lambdas(3)+lambdas(4))*M, obj.Pfem);
       %     obj.qnodal=obj.fromFEMVector(obj.qfem(:,1));
       %     qfem=obj.qfem;
       % end

       function conds = tabMatrixCondition( obj, alphas, mode )
           conds=[];
           disp(['Matrix conditions for mode: ' num2str(mode)]);
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           K  =  obj.assemblyGlobalMatrix('computeStifnessMatrix', 1, obj.isConst);
           M  =  obj.assemblyGlobalMatrix('computeMassMatrix', 1, obj.isConst);
           [~, lambdas] = solver.solveEigenproblem(K,M,10);
           lambdas = diag(lambdas(1:5,1:5));

           for k=1:numel(alphas)
               disp(['iter: ' num2str(k)]);
               conds = [conds rcond(full(solver.createSparseMatrix(K-alphas(k)*lambdas(mode)*M)))];
           end
       end
       
   end
end


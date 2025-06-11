classdef ElasticHarmonicVibrations < FEAnalysis
   
   properties
        isMeshConst, isLoadConst, lambda, omegas, mode, modes, count, P0;
   end
   
   methods       
       
       function obj = ElasticHarmonicVibrations(felems, mesh, mode, isMeshConst, isLoadConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.mode = mode;
            obj.isMeshConst=isMeshConst;
            obj.isLoadConst=isLoadConst;
            obj.rotations=[];
            obj.count=0;
            obj.modes=[];
       end

       function [K, M] = computeMatrices(obj,x)
           K  = obj.createSparseMatrix( obj.assemblyGlobalMatrix('computeStifnessMatrix', x, obj.isMeshConst) );
           M  = obj.createSparseMatrix( obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isMeshConst) );
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

       %Harmonic Vibrations
       function qfem = solve(obj, x)
           obj.count=obj.count+1;
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           K  = obj.assemblyGlobalMatrix('computeStifnessMatrix', x, obj.isMeshConst);
           M  = obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isMeshConst);

           [modes, lambdas] = solver.solveEigenproblem(K,M,10);
           obj.modes = cat(3, obj.modes, modes);
           lambdas = diag(lambdas(1:5,1:5));
           obj.omegas=sqrt(lambdas/2/pi);
           obj.Pfem=0*obj.Pfem;
           %obj.qfem = solver.solve(K-(lambdas(3)+lambdas(4))/2*M, obj.Pfem);
           obj.Pfem = obj.toFEMVector(lambdas(obj.mode)*obj.fromFEMVector(modes(:,obj.mode)));
           if obj.count==1
               obj.P0=obj.Pfem;
           end
           if obj.isLoadConst
                obj.qfem = solver.solve(K, obj.P0);
           else
                obj.qfem = solver.solve(K, obj.Pfem);
           end
           obj.qnodal=obj.fromFEMVector(obj.qfem(:,1));
           qfem=obj.qfem;
       end

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
        
       % function corr=computeCorrelations(nmodes)
       %     corr=zeros(1,nmodes);
       %     for k=1:nmodes
       %          corr(k)=abs(mode(:,k-1)'*mode(:,k))/norm(mode(:,k-1))/norm(mode(:,k))
       %     end
       % end
       % 
   end
end


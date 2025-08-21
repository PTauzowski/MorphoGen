classdef ElasticHarmonicVibrations < FEAnalysis
   
   properties
        isMeshConst, isLoadConst, lambdas, omegas, frequencies, mode, modes, nmodes, count, P0, correlation_matrix, freedofs;
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
            obj.lambdas=[];
            obj.omegas=[];
            obj.frequencies=[];
            obj.nmodes=30;
            obj.correlation_matrix=zeros(obj.nmodes,obj.nmodes);
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
           vK  = obj.assemblyGlobalMatrix('computeStifnessMatrix', x, obj.isMeshConst);
           vM  = obj.assemblyGlobalMatrix('computeMassMatrix', x, obj.isMeshConst);
           obj.freedofs=solver.freedofs;

           K = solver.createSparseMatrix(vK);
           M = solver.createSparseMatrix(vM);
           save("LoadForce.mat",  "K", "M", "x");

           [modes, lambdas] = solver.solveEigenproblem(vK,vM,obj.nmodes);
           obj.modes = cat(3, obj.modes, modes);
           obj.lambdas = [obj.lambdas diag(lambdas)];
           obj.omegas=[ obj.omegas sqrt(diag(lambdas))];
           obj.frequencies = [ obj.frequencies sqrt(diag(lambdas)) / 2 / pi];
           obj.Pfem=0*obj.Pfem;
           %obj.qfem = solver.solve(K-(lambdas(3)+lambdas(4))/2*M, obj.Pfem);
           Pfem=reshape(modes(solver.freedofs,obj.mode)),size(obj.ndofs,2),size(solver.freedofs,1))';
           obj.Pfem = obj.toFEMVector(lambdas(obj.mode,obj.mode)*M*obj.fromFEMVector(modes(solver.freedofs,obj.mode)));
           if obj.count==1
               obj.P0=obj.Pfem;
           end
           if obj.isLoadConst
                obj.qfem = solver.solve(vK, obj.P0);
           else
                obj.qfem = solver.solve(vK, obj.Pfem);
           end
           obj.qnodal=obj.fromFEMVector(obj.qfem(:,1));
           qfem=obj.qfem;
       end

       function Pfem = computeLoadVector(obj, frame)
           Pfem = obj.toFEMVector(obj.lambdas(obj.mode,frame)*M*obj.fromFEMVector(obj.modes(obj.freedofs,obj.mode,frame)));
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
        
       function computeSelfCorrelationMatrix(obj)
           for k=1:obj.nmodes
               for l=1:obj.nmodes
                    obj.correlation_matrix(k,l)=abs(obj.modes(:,k)'*obj.modes(:,l))/norm(obj.modes(:,k))/norm(obj.modes(:,l));
               end
           end
       end

       function Corr = computeCorrelationMatrix(obj, dd_modes, nmodes)
           Corr=zeros(nmodes,nmodes);
           for k=1:nmodes
               for l=1:nmodes
                    Corr(k,l)=abs(dd_modes(:,k)'*obj.modes(:,l,end))/norm(dd_modes(:,k))/norm(obj.modes(:,l,end));
               end
           end
       end

       function cv = computeCorrelationVector(obj,test_mode)
           cv=zeros(1,obj.nmodes);
           for k=1:obj.nmodes
                    cv(k)=abs(test_mode'*obj.modes(:,k,end))/norm(test_mode)/norm(obj.modes(:,k,end));
           end
       end

   end
end


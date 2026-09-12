classdef SecondOrderDynamicElasticityWeighted < FEAnalysis
   
   properties
        isConst, lambda, omegas, lambdac, modes;
   end
   
   methods       
       function obj = SecondOrderDynamicElasticityWeighted(felems, mesh, lambdac, isConst)
            obj= obj@FEAnalysis( felems, mesh );
            obj.lambdac = lambdac;
            obj.isConst=isConst;
            obj.rotations=[];
            obj.modes=[];

       end
       
       function qfem = solve(obj, x)
           [I,J,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
            if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
           end
           if obj.isConst
               K = obj.assemblyGlobalMatrix('computeStifnessMatrixConst',x,false);
           else
               K = obj.assemblyGlobalMatrix('computeStifnessMatrix',x,false);
           end
           obj.qfem = solver.solve(K, obj.Pfem);
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults(x);
           Kg = obj.assemblyGlobalMatrix('computeGeometricStifnessMatrix',x,false);
           M  = obj.assemblyGlobalMatrix('computeMassMatrix',x,false);
           [eigenvectors, lambdas] = solver.solveEigenproblem(K,Kg,10);
           [mode, omegas1] = solver.solveEigenproblem(K,M,10);
           omegas = sqrt(omegas1);
           obj.modes = cat(3, obj.modes, mode);  
           obj.lambda=lambdas(1);
           obj.omegas=diag(omegas(1:5,1:5));
           obj.qfem=0*obj.Pfem;
           %obj.qfem(solver.freedofs)=eigenvectors(:,1);
           obj.qfem = solver.solve(K-abs(omegas(1))*Kg, obj.Pfem);
           obj.qnodal=obj.fromFEMVector(obj.qfem(:,1));
           qfem=obj.qfem;

       end

       
       
   end
end

